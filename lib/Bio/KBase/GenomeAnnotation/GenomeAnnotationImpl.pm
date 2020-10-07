package Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use strict;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

GenomeAnnotation

=head1 DESCRIPTION

API Access to the Genome Annotation Service.

  Provides support for gene calling, functional annotation, re-annotation. Use to extract annotation in
formation about an existing genome, or to create new annotations.

=cut

#BEGIN_HEADER

use Safe;
use File::Basename;
use File::Temp;
use File::Which;
use File::Slurp;
use Data::Dumper;
use Digest::MD5 'md5_hex';
use Time::HiRes 'gettimeofday';
use POSIX;
use SOAP::Lite;
use Clone;

our $have_kbase_idserver;
eval {
    require Bio::KBase::IDServer::Client;
    $have_kbase_idserver = 1;
};

# disable # use Bio::KBase::KmerAnnotationByFigfam::Client;
# disable # use KmerClassifier;
#use Bio::KBase::KIDL::Helpers qw(json_to_tempfile tempfile_to_json);
use IPC::Run qw(run);
use JSON::XS;
use File::Slurp;
use Bio::KBase::GenomeAnnotation::Awe;
use Bio::KBase::GenomeAnnotation::Shock;

use Bio::KBase::GenomeAnnotation::Glimmer;
use BinningReports;
use GenomeTypeObject;
use gjogenbank;
use GenBankToGTO;
use IDclient;
# disabled # use ANNOserver;
use SeedUtils;
use SeedURLs;
use gjoseqlib;
use StrepRepeats;
use overlap_resolution;
use PropagateGBMetadata;
use Capture::Tiny 'capture_stderr';
use AdaboostClassify;
use MongoDB;

use Bio::KBase::GenomeAnnotation::SubsystemProjector;
use Bio::P3::DeploymentConfig;

sub _get_coder
{
    return JSON::XS->new->ascii->pretty->allow_nonref;
}

sub _call_using_strep_repeats
{
    my($self, $genome_in, $tool) = @_;
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $parsed = StrepRepeats::get_raw_repeats($genome_in, $tool);

    my $event = {
	tool_name => $tool,
	execution_time => scalar gettimeofday,
	parameters => [],
	hostname => $self->{hostname},
    };

    my $idc = IDclient->new($genome_in);
    my $event_id = $genome_in->add_analysis_event($event);

    my $type = 'repeat_unit';

    for my $r (@$parsed)
    {
	my $loc = $r->{location};
	$genome_in->add_feature({
	    -id_client 	     => $idc,
	    -id_prefix 	     => $genome_in->{id},
	    -type 	     => $type,
	    -location 	     => $loc,
	    -function 	     => $r->{function},
	    -annotation      => $r->{annotation},
	    -analysis_event_id 	     => $event_id,
	});
    }

    $genome_in = $genome_in->prepare_for_return();

    return $genome_in;
}

sub _allocate_local_genome_id
{
    my($self, $taxon_id, $mongo_host, $mongo_db) = @_;

    my $query_timeout = -1;
    my $connection_timeout = 120;

    $mongo_db //= "seed-genome-allocation";
 
    my $client = MongoDB->connect($mongo_host);
    #my $client = MongoDB::MongoClient->new(host => $mongo_host, query_timeout => $query_timeout,
#					   timeout => $connection_timeout,
#					   auto_reconnect => 1,
#					   auto_connect => 1);
#
    my $db = $client->get_database($mongo_db);
    my $coll_next = $db->get_collection('next');

    $coll_next->ensure_index({prefix => 1});

    my $res = $coll_next->find_one_and_update({ taxon => $taxon_id }, { '$inc' => { next_val => 1 } }, { upsert => 1, new => 1 });
    #print Dumper($res);

    if (!ref($res))
    {
	die "MongoDB error: $res";
    }

    my $val = $res->{next_val};

    return "$taxon_id.0$val";
}


sub _allocate_seed_genome_id
{
    my($self, $taxon_id, $url) = @_;

    my $proxy = SOAP::Lite->uri('http://www.soaplite.com/Scripts')-> proxy($url);

    #
    # If we have heavy registration traffic we may fail with messages like
    #
    #   Error registering genome via SEED clearinghouse: soap:Server error executing queries: DBD::mysql::db do failed: Deadlock found when trying to get lock; try restarting transaction at /vol/core-seed/Clearinghouse/FIG/CGI/clearinghouse_services.cgi line 162.
    #
    # So we will do a backoff and retry here.
    #

    my $max_retries = 10;
    my $try = 0;
    while ($try++ < $max_retries)
    {
	my $r = $proxy->register_genome($taxon_id);
	if ($r->fault) {
	    warn "Error on try $try registering genome via SEED clearinghouse: " . $r->faultcode . " " . $r->faultstring;
	    sleep(2 + rand(10));
	}
	else
	{
	    my $id = $r->result;
	    return "$taxon_id.$id";
	}
    }
    #
    # Retries failed
    #
    die "Unable to register genome via SEED clearinghouse after $max_retries tries\n";
}

sub _allocate_kb_genome_id
{
    my($self, $taxon_id, $url) = @_;

    my $id_prefix = "kb";
	
    my $idc = Bio::KBase::IDServer::Client->new($url);
    my $id = "$id_prefix|g." . $idc->allocate_id_range("$id_prefix|g", 1);
    return $id;
}

sub _allocate_genome_id_fallback
{
    my($self, $taxon_id) = @_;

    my $id = "random|$taxon_id." . int(rand(10000));
    return $id;
}

#
# Lazy load of projection data.
#
# This may return undef if the projector was disabled at service startup.
#
# If it gets too expensive, support loading of serialized objet via e.g. Sereal.
#
sub _subsystem_projector
{
    my($self) = @_;
    if (exists $self->{_subsystem_projector})
    {
	return $self->{_subsystem_projector};
    }

    my $proj = Bio::KBase::GenomeAnnotation::SubsystemProjector->new($self->{subsystem_roles},
								     $self->{subsystem_variants});
    $proj->load_reference_subsystems($self->{subsystem_reference_data});
    $proj->load_variant_codes($self->{subsystem_variant_map});
    $self->{_subsystem_projector} = $proj;
    return $proj;
}

#
# Some gene callers (phanotate) may return coordinates of the form "<N". For now
# patch these by replacing with N-1 (if N>1) or 1.
#
sub _fix_fuzzy_coordinate
{
    my($n, $contig_length) = @_;

    if ($n =~ /<(\d+)/)
    {
	$n = $1 - 1;
	$n = 1 if $n < 1;
    }
    elsif ($n =~ />(\d+)/)
    {
	$n = $1 + 1;
	$n = $contig_length if $contig_length && $n > $contig_length;
    }

    return $n;
}

#
# Allocate a new feature ID in the given genome with an initialized id client.
#
sub _allocate_new_feature_id
{
    my($idc, $from_id, $type) = @_;
    my $prefix = $from_id;
    $prefix =~ s/\.\d+$//;
    my $suffix = $idc->allocate_id_range($prefix, 1);
    my $id = join(".", $prefix, $suffix);
    return $id;
}
#
# translate strand/left/right to location
#
sub _ends_to_location
{
    my($ctg, $strand, $left, $right) = @_;
    if ($strand eq '+')
    {
	return [$ctg, $left, $strand, $right - $left + 1];
    }
    else
    {
	return [$ctg, $right, $strand, $right - $left + 1];
    }
}

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $cfg = Bio::P3::DeploymentConfig->new($ENV{KB_SERVICE_NAME} || "GenomeAnnotation");

    my $dir = $cfg->setting("kmer_v2_data_directory");
    $dir or die "Configuration parameter for kmer_v2_data_directory not set";
    -d $dir or die "Directory $dir for kmer_v2_data_directory does not exist";
	
    $self->{kmer_v2_data_directory} = $dir;

    if (my $temp = $cfg->setting("tempdir"))
    {
	$ENV{TEMPDIR} = $ENV{TMPDIR} = $temp;
    }

    my $dir = $cfg->setting("kmer_classifier_data_directory");
    #
    # Make these soft errors for now.
    # 2018-0827 - silence them since this is functionality that never really came online
    #
    if (0)
    {
	if (!$dir)
	{
	    warn "Configuration parameter for kmer_classifier_data_directory not set";
	}
	elsif (! -d $dir)
	{
	    warn "Directory $dir for kmer_classifier_data_directory does not exist";
	}
    }
	
    $self->{kmer_classifier_data_directory} = $dir;

    $dir = $cfg->setting("nr-annotation-directory");
    if (!$dir)
    {
	warn "Configuration parameter for nr-annotation-directory not set";
    }
    elsif (! -d $dir)
    {
	warn "Directory $dir for nr-annotation_directory does not exist";
    }
	
    $self->{nr_annotation_directory} = $dir;

    $self->{vigor_reference_db_directory} = $cfg->setting("vigor-reference-db-directory");

    #
    # This is stand-in code; read the viral family data into memory.
    #
    my $vfam = $cfg->setting("viral-family-db");
    $self->{viral_plfam} = {};
    $self->{viral_pgfam} = {};
    if ($vfam && open(my $vf, "<", $vfam))
    {
	$_ = <$vf>;
	while (<$vf>)
	{
	    chomp;
	    my($type, $id, $product) = split(/\t/);
	    if ($type eq 'plfam')
	    {
		my($genus) = $id =~ /^PLF_(\d+)_/;
		$self->{viral_plfam}->{$genus}->{$product} = [$type, $id, $product];
	    }
	    elsif ($type eq 'pgfam')
	    {
		$self->{viral_pgfam}->{$product} = [$type, $id, $product];
	    }
	}
	close($vf);
    }

    my $phage_files = $cfg->setting("phage-annotation-files");
    if (!$phage_files)
    {
	$phage_files = [];
    }
    elsif (!ref($phage_files))
    {
	$phage_files = [$phage_files];
    }
    $self->{phage_annotation_files} = $phage_files;

    $self->{special_protein_threads} = $cfg->setting("special_protein_threads") // 1;
    $self->{special_protein_dbdir} = $cfg->setting("special_protein_dbdir");
    $self->{special_protein_lookup_db} = $cfg->setting("special_protein_lookup_db");
    $self->{special_protein_cache_db} = $cfg->setting("special_protein_cache_db");
    $self->{special_protein_cache_dbhost} = $cfg->setting("special_protein_cache_dbhost");
    $self->{special_protein_cache_dbuser} = $cfg->setting("special_protein_cache_dbuser");
    $self->{special_protein_cache_dbpass} = $cfg->setting("special_protein_cache_dbpass");

    $self->{patric_call_proteins_remote_host} = $cfg->setting("patric_call_proteins_remote_host");
    $self->{patric_call_proteins_remote_user} = $cfg->setting("patric_call_proteins_remote_user");
    $self->{patric_call_proteins_remote_key} = $cfg->setting("patric_call_proteins_remote_key");
    $self->{patric_call_proteins_path} = $cfg->setting("patric_call_proteins_path");
    $self->{patric_call_proteins_ff_path} = $cfg->setting("patric_call_proteins_ff_path");
    $self->{patric_call_proteins_md5_to_fam_path} = $cfg->setting("patric_call_proteins_md5_to_fam_path");

    $self->{patric_annotate_families_url} = $cfg->setting("patric_annotate_families_url");
    $self->{patric_annotate_families_kmers} = $cfg->setting("patric_annotate_families_kmers");

    $self->{mongo_allocate_id_host} = $cfg->setting("mongo-allocate-id-host");
    $self->{mongo_allocate_id_db} = $cfg->setting("mongo-allocate-id-db");

    $self->{patric_mlst_dbdir} = $cfg->setting("patric_mlst_dbdir");
    $self->{workflow_dir} = $cfg->setting("workflow-dir");

    #
    # Evaluation data files.
    #

    $self->{genome_evaluation_data} = $cfg->setting('genome-evaluation-data');
    $self->{genome_evaluation_predictors} = $cfg->setting('genome-evaluation-predictors');
    $self->{genome_evaluation_checkg} = $cfg->setting('genome-evaluation-checkg');
    $self->{seedtk_path} = $cfg->setting('seedtk-path');

    #
    # String to use in CDS identifiers.
    #
    $self->{cds_id_type} = $cfg->setting("cds_id_type") || "CDS";

    #
    # Determine if we are using KB or SEED genome ID generation.
    #

    my $surl = $cfg->setting("seed-id-clearinghouse");
    my $iurl = $cfg->setting("idserver_url");

    #
    # For backward compatibility set this
    #
    $iurl = 'https://kbase.us/services/idserver' if !defined($iurl);

    if ($self->{mongo_allocate_id_host})
    {
	print STDERR "Using local mongo for id allocation\n";
	$self->{allocate_genome_id} = sub {
	    my($taxon_id) = @_;
	    $self->_allocate_local_genome_id($taxon_id, $self->{mongo_allocate_id_host}, $self->{mongo_allocate_id_db});
	};
    }
    elsif ($surl)
    {
	print STDERR "Using SEED id allocation from $surl\n";
	$self->{allocate_genome_id} = sub { my($taxon_id) = @_;
					    $self->_allocate_seed_genome_id($taxon_id, $surl);
					};
    }
    elsif ($iurl)
    {
	$have_kbase_idserver or die "Using KBase id allocation but KBase ID server code not available";
	print STDERR "Using KBase id allocation from $iurl\n";
	$self->{allocate_genome_id} = sub { my($taxon_id) = @_;
					    $self->_allocate_kb_genome_id($taxon_id, $iurl);
					};
    }
    else
    {
	print STDERR "No id allocation mechanism defined\n";
	$self->{allocate_genome_id} = sub { my($taxon_id) = @_;
					    $self->_allocate_genome_id_fallback($taxon_id);
					};
    }
    
    $self->{kmer_service_url} = $cfg->setting("kmer_service_url");
    $self->{awe_server} = $cfg->setting("awe-server");
    $self->{shock_server} = $cfg->setting("shock-server");
    $self->{genemark_home} = $cfg->setting("genemark-home");

    print STDERR "kmer_v2_data_directory = $self->{kmer_v2_data_directory}\n";

    #
    # Subsystem propagation data files.
    # If any are defined and missing on disk, croak.
    # If any are not defined, croak.
    # If all are defined and present, proceed.
    # If all are missing, proceed, and flag subsystem propagation as not present.
    #
    {
	my @keys = qw(subsystem-roles subsystem-variants subsystem-reference-data subsystem-variant-map);
	my %settings;
	my @missing;
	my @not_found;
	for my $k (@keys)
	{
	    my $f = $cfg->setting($k);
	    print "Subsystem setting '$k': '$f'\n";
	    if ($f)
	    {
		if (-f $f)
		{
		    my $sk = $k;
		    $sk =~ s/-/_/g;
		    $settings{$sk} = $f;
		}
		else
		{
		    push @not_found, [$k, $f];
		}
	    }
	    else
	    {
		push(@missing, $k);
	    }
	}
	if (scalar @keys == scalar @missing)
	{
	    print STDERR "Subsystem projection data not defined, disabling.\n";
	    $self->{_subsystem_projector} = undef;
	}
	elsif (@not_found || @missing)
	{
	    die "Subsystem projection data files not found:\n" . join("\n", map { "\t$_->[0]: $_->[1]" } @not_found) . "\nOr missing: @missing";
	}
	$self->{$_} = $settings{$_} foreach keys %settings;
    }

    my $h = `hostname`;
    chomp $h;
    $self->{hostname} = $h;
    
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}
=head1 METHODS
=head2 genome_ids_to_genomes

  $genomes = $obj->genome_ids_to_genomes($ids)

=over 4




=item Description

Given one or more Central Store genome IDs, convert them into genome objects.
=back

=cut

sub genome_ids_to_genomes
{
    my $self = shift;
    my($ids) = @_;

    my @_bad_arguments;
    (ref($ids) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"ids\" (value was \"$ids\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to genome_ids_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genomes);
    #BEGIN genome_ids_to_genomes

    

    #END genome_ids_to_genomes
    my @_bad_returns;
    (ref($genomes) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"genomes\" (value was \"$genomes\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to genome_ids_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genomes);
}


=head2 create_genome

  $genome = $obj->create_genome($metadata)

=over 4




=item Description

Create a new genome object and assign metadata.
=back

=cut

sub create_genome
{
    my $self = shift;
    my($metadata) = @_;

    my @_bad_arguments;
    (ref($metadata) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"metadata\" (value was \"$metadata\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome

    $genome = GenomeTypeObject->new();
    $genome->set_metadata($metadata);

    #
    # If no ID was created, allocate a genome ID from our ID server.
    #
    if (!$genome->{id})
    {
	my $tax = $genome->{ncbi_taxonomy_id};
	
	if ($tax !~ /^\d+$/)
	{
	    $tax = "6666666";
	}

	$genome->{id} = $self->{allocate_genome_id}->($tax);
    }
    $genome = $genome->prepare_for_return();

    #END create_genome
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome);
}


=head2 create_genome_from_genbank

  $genome = $obj->create_genome_from_genbank($gb_data)

=over 4




=item Description

Create a new genome object from one or more genbank files.
=back

=cut

sub create_genome_from_genbank
{
    my $self = shift;
    my($gb_data) = @_;

    my @_bad_arguments;
    (!ref($gb_data)) or push(@_bad_arguments, "Invalid type for argument \"gb_data\" (value was \"$gb_data\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome_from_genbank:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome_from_genbank

    #
    # We need to parse the file into a list of entries first so that
    # we may extract the taxon ID and allocate a genome ID. This way
    # our initial feature set can be created with the appropriate identifiers.
    #

    my @entries = parse_genbank(\$gb_data);
    
    my $tax_id;
    for my $entry (@entries)
    {
	my @sources = gjogenbank::features_of_type( $entry, 'source' );
	my ( $src_taxid ) = map { /^taxon:(\S+)$/ ? $1 : () }
		map { $_->[1]->{db_xref} ? @{$_->[1]->{db_xref}} : () }
		@sources;
	if ($src_taxid)
	{
	    $tax_id = $src_taxid;
	    last;
	}
    }

    if ($tax_id !~ /^\d+$/)
    {
	$tax_id = "6666666";
    }
    
    my $genome_id = $self->{allocate_genome_id}->($tax_id);

    $genome = GenBankToGTO::new({ entry => \@entries, id => $genome_id });

    #
    # The current genbank code is creating features with type peg; remap them to CDS.
    #

    for my $feature ($genome->features)
    {
	if ($feature->{type} eq 'peg')
	{
	    $feature->{type} = 'CDS';
	}
    }

    $genome->prepare_for_return();

    #END create_genome_from_genbank
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome_from_genbank:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome);
}


=head2 create_genome_from_SEED

  $genome = $obj->create_genome_from_SEED($genome_id)

=over 4




=item Description

Create a new genome object based on data from the SEED project.
=back

=cut

sub create_genome_from_SEED
{
    my $self = shift;
    my($genome_id) = @_;

    my @_bad_arguments;
    (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_id\" (value was \"$genome_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome_from_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome_from_SEED
    #END create_genome_from_SEED
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome_from_SEED:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome);
}


=head2 create_genome_from_RAST

  $genome = $obj->create_genome_from_RAST($genome_or_job_id)

=over 4




=item Description

Create a new genome object based on a RAST genome.
=back

=cut

sub create_genome_from_RAST
{
    my $self = shift;
    my($genome_or_job_id) = @_;

    my @_bad_arguments;
    (!ref($genome_or_job_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_or_job_id\" (value was \"$genome_or_job_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome_from_RAST:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome_from_RAST

    print STDERR "get RAST : ctx=" . Dumper($ctx);
    $genome = {};
    #END create_genome_from_RAST
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome_from_RAST:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome);
}


=head2 set_metadata

  $genome_out = $obj->set_metadata($genome_in, $metadata)

=over 4




=item Description

Modify genome metadata.
=back

=cut

sub set_metadata
{
    my $self = shift;
    my($genome_in, $metadata) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($metadata) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"metadata\" (value was \"$metadata\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN set_metadata

    $genome_out = GenomeTypeObject->initialize_without_indexes($genome_in);
    $genome_out->set_metadata($metadata);
    $genome_out = $genome_out->prepare_for_return();
    
    #end set_metadata
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'set_metadata');
    }
    return($genome_out);
}




=head2 add_contigs

  $genome_out = $obj->add_contigs($genome_in, $contigs)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Add a set of contigs to the genome object.

=back

=cut

sub add_contigs
{
    my $self = shift;
    my($genome_in, $contigs) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"contigs\" (value was \"$contigs\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_contigs');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #END set_metadata
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 add_contigs

  $genome_out = $obj->add_contigs($genome_in, $contigs)

=over 4




=item Description

Add a set of contigs to the genome object.
=back

=cut

sub add_contigs
{
    my $self = shift;
    my($genome_in, $contigs) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"contigs\" (value was \"$contigs\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN add_contigs

    $genome_out = GenomeTypeObject->initialize_without_indexes($genome_in);
    $genome_out->add_contigs($contigs);
    $genome_out = $genome_out->prepare_for_return();

    #END add_contigs
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 add_contigs_from_handle

  $genome_out = $obj->add_contigs_from_handle($genome_in, $contigs)

=over 4




=item Description

Add a set of contigs to the genome object, loading the contigs
from the given handle service handle.
=back

=cut

sub add_contigs_from_handle
{
    my $self = shift;
    my($genome_in, $contigs) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"contigs\" (value was \"$contigs\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_contigs_from_handle:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN add_contigs_from_handle
    #END add_contigs_from_handle
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to add_contigs_from_handle:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 add_features

  $genome_out = $obj->add_features($genome_in, $features)

=over 4




=item Description

Add a set of features in tabular form.
=back

=cut

sub add_features
{
    my $self = shift;
    my($genome_in, $features) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($features) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"features\" (value was \"$features\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN add_features

    $genome_out = GenomeTypeObject->initialize($genome_in);
    $genome_out->add_features_from_list($features);
    $genome_out = $genome_out->prepare_for_return();
    
    #END add_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to add_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 genomeTO_to_reconstructionTO

  $return = $obj->genomeTO_to_reconstructionTO($genomeTO)

=over 4




=item Description


=back

=cut

sub genomeTO_to_reconstructionTO
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to genomeTO_to_reconstructionTO:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN genomeTO_to_reconstructionTO

    die "method genomeTO_to_reconstructionTO is disabled";
    use ANNOserver;
    my $annoO = ANNOserver->new;
    my %in_models = map { chop; ($_ => 1) } `all_roles_used_in_models`;
    my %bindings_to_roles;

    my $features = $genomeTO->{features};
    my @role_fid_tuples;
    my $assignments = [];
    foreach my $fidH (@$features)
    {
	my $fid = $fidH->{id};
	my $f = $fidH->{function};
	if ($f)
	{
	    
	    push(@$assignments,[$fid,$f]);
	    foreach my $role (&SeedUtils::roles_of_function($f))
	    {
		push(@role_fid_tuples,[$role,$fid]);
		if ($in_models{$role}) { $bindings_to_roles{$role}->{$fid} = 1 }
	    }
	}
    }
    my $mr = $annoO->metabolic_reconstruction({-roles => \@role_fid_tuples});
    my $sub_vars = [];
    my $bindings = [];
    my %subsys;
    foreach my $tuple (@$mr)
    {
	my($sub_var,$role,$fid) = @$tuple;
	my($sub,$var) = split(/:/,$sub_var);
	if ($var !~ /\*?(0|-1)\b/)
	{
	    $subsys{$sub} = $var;
	    $bindings_to_roles{$role}->{$fid} = 1;
	}
    }
    foreach my $role (keys(%bindings_to_roles))
    {
	my $roles = $bindings_to_roles{$role};
	my @fids = keys(%$roles);
	foreach my $fid (@fids)
	{
	    push(@$bindings,[$fid,$role]);
	}
    }
    my @sv = map { [$_,$subsys{$_}] } keys(%subsys);
    $return = {
	subsystems => \@sv,
	assignments => $assignments,
	bindings   => $bindings,
    };
    
    #END genomeTO_to_reconstructionTO
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to genomeTO_to_reconstructionTO:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 genomeTO_to_feature_data

  $return = $obj->genomeTO_to_feature_data($genomeTO)

=over 4




=item Description


=back

=cut

sub genomeTO_to_feature_data
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to genomeTO_to_feature_data:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN genomeTO_to_feature_data

    my $feature_data = [];
    my $features = $genomeTO->{features};
    foreach my $feature (@$features)
    {
	my $fid = $feature->{id};
	my $loc = join(",",map { my($contig,$beg,$strand,$len) = @$_; 
				 "$contig\_$beg$strand$len" 
			       } @{$feature->{location}});
	my $type = $feature->{type};
	my $func = $feature->{function};
	my $md5 = "";
	$md5 = md5_hex(uc($feature->{protein_translation})) if $feature->{protein_translation};
	my $aliases = join(",",@{$feature->{aliases} // []});
	push(@$feature_data,[$fid,$loc,$type,$func,$aliases,$md5]);
    }
    $return = $feature_data;
    #END genomeTO_to_feature_data
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to genomeTO_to_feature_data:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 reconstructionTO_to_roles

  $return = $obj->reconstructionTO_to_roles($reconstructionTO)

=over 4




=item Description


=back

=cut

sub reconstructionTO_to_roles
{
    my $self = shift;
    my($reconstructionTO) = @_;

    my @_bad_arguments;
    (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"reconstructionTO\" (value was \"$reconstructionTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to reconstructionTO_to_roles:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN reconstructionTO_to_roles

    my $bindings = $reconstructionTO->{bindings};
    my %roles = map { ($_->[1] => 1) } @$bindings;
    $return = [sort keys(%roles)];

    #END reconstructionTO_to_roles
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to reconstructionTO_to_roles:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 reconstructionTO_to_subsystems

  $return = $obj->reconstructionTO_to_subsystems($reconstructionTO)

=over 4




=item Description


=back

=cut

sub reconstructionTO_to_subsystems
{
    my $self = shift;
    my($reconstructionTO) = @_;

    my @_bad_arguments;
    (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"reconstructionTO\" (value was \"$reconstructionTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to reconstructionTO_to_subsystems:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN reconstructionTO_to_subsystems

    my $subsys_pairs = $reconstructionTO->{subsystems};
    $return = $subsys_pairs;
    
    #END reconstructionTO_to_subsystems
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to reconstructionTO_to_subsystems:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 assign_functions_to_CDSs

  $return = $obj->assign_functions_to_CDSs($genomeTO)

=over 4




=item Description

Given a genome object populated with contig data, perform gene calling
and functional annotation and return the annotated genome.
=back

=cut

sub assign_functions_to_CDSs
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to assign_functions_to_CDSs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN assign_functions_to_CDSs

    die "method assign_functions_to_CDSs is disabled";
    #END assign_functions_to_CDSs
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to assign_functions_to_CDSs:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 annotate_genome

  $return = $obj->annotate_genome($genomeTO)

=over 4




=item Description


=back

=cut

sub annotate_genome
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_genome

    $return = $self->run_pipeline($genomeTO, $self->default_workflow());
    
    #END annotate_genome
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_selenoproteins

  $return = $obj->call_selenoproteins($genomeTO)

=over 4




=item Description


=back

=cut

sub call_selenoproteins
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_selenoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_selenoproteins

    my $idc = IDclient->new($genomeTO);

    my $coder = _get_coder();
    
    my $genomeTO_json = $coder->encode($genomeTO);

    my $genomeOut_json;

    my $tmp = File::Temp->new();
    print $tmp $genomeTO_json;
    close($tmp);

    my $id_prefix = $genomeTO->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }

    my @cmd = ('rast_call_special_proteins',
	       '--seleno',
	       "--id-type", $self->{cds_id_type},
	       '--id-prefix', $id_prefix,
	       '--input', $tmp);
    $ctx->stderr->log_cmd(@cmd);
    my $ok = run(\@cmd,
		 '>', \$genomeOut_json,
		 $ctx->stderr->redirect);

    undef $tmp;
    undef $genomeTO;

    if ($ok)
    {
	$return = $coder->decode($genomeOut_json);
    }
    else
    {
	die "rast_call_special_proteins failed: ", $ctx->stderr->text_value;
    }
    
    #END call_selenoproteins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_selenoproteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_pyrrolysoproteins

  $return = $obj->call_pyrrolysoproteins($genomeTO)

=over 4




=item Description


=back

=cut

sub call_pyrrolysoproteins
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_pyrrolysoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_pyrrolysoproteins
    
    my $coder = _get_coder();
    
    my $genomeTO_json = $coder->encode($genomeTO);

    my $genomeOut_json;

    my $tmp = File::Temp->new();
    print $tmp $genomeTO_json;
    close($tmp);

    my $id_prefix = $genomeTO->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }

    my @cmd = ('rast_call_special_proteins',
	       '--pyrro',
	       "--id-type", $self->{cds_id_type},
	       '--id-prefix', $id_prefix,
	       '--input', $tmp);
    $ctx->stderr->log_cmd(@cmd);
			  
    my $ok = run(\@cmd,
		 '>', \$genomeOut_json,
		 $ctx->stderr->redirect);

    undef $tmp;
    undef $genomeTO;

    if ($ok)
    {
	$return = $coder->decode($genomeOut_json);
    }
    else
    {
	die "rast_call_special_proteins failed: " . $ctx->stderr->text_value;
    }

    #END call_pyrrolysoproteins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_pyrrolysoproteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_selenoprotein

  $return = $obj->call_features_selenoprotein($genomeTO)

=over 4




=item Description

Given a genome typed object, call selenoprotein features.
=back

=cut

sub call_features_selenoprotein
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_selenoprotein:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_selenoprotein

    $return = $self->call_selenoproteins($genomeTO);
    
    #END call_features_selenoprotein
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_selenoprotein:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_pyrrolysoprotein

  $return = $obj->call_features_pyrrolysoprotein($genomeTO)

=over 4




=item Description

Given a genome typed object, call pyrrolysoprotein features.
=back

=cut

sub call_features_pyrrolysoprotein
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_pyrrolysoprotein:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_pyrrolysoprotein

    $return = $self->call_pyrrolysoproteins($genomeTO);
    
    #END call_features_pyrrolysoprotein
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_pyrrolysoprotein:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_insertion_sequences

  $return = $obj->call_features_insertion_sequences($genomeTO)

=over 4




=item Description

Given a genome typed object, call insertion sequences.
=back

=cut

sub call_features_insertion_sequences
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_insertion_sequences:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_insertion_sequences

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
    my $output_file = File::Temp->new();

    my @params = ();
    my @cmd = ("find_IS", @params);
    $ctx->stderr->log_cmd(@cmd);

    my $ok = run(\@cmd,
		 "<", $sequences_file,
		 ">", $output_file,
		 $ctx->stderr->redirect);

    unlink($sequences_file);

    if (!$ok)
    {
	die "Error running find_IS: @cmd\n";
    }

    close($output_file);
    my($res_fh);
    open($res_fh, "<", $output_file) or die "Cannot open find_IS output file $output_file: $!";

    my $event = {
	tool_name => $cmd[0],
	execution_time => scalar gettimeofday,
	parameters => \@params,
	hostname => $self->{hostname},
    };

    my $idc = IDclient->new($genome_in);
    my $event_id = $genome_in->add_analysis_event($event);

    #	fig|NC_002952.IS.1	NC_002952_760426_758423	IS1182	IS1182 family insertion sequence
    #	fig|NC_002952.IS.2	NC_002952_2545636_2543556	IS_Unknown	Putative insertion sequence

    my $id_prefix = $genome_in->{id};

    while(<$res_fh>)
    {
	chomp;
	my($id, $loc, $type, $function) = split(/\t/);

	my $typed_prefix = join(".", $id_prefix, $type);

	my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, 1);

	my $id = join(".", $typed_prefix, $cur_id_suffix);

	my $locs = GenomeTypeObject::seed_location_to_location_list($loc);

	$genome_in->add_feature({
	    -id              => $id,
	    -type 	     => $type,
	    -location 	     => $locs,
	    -function 	     => $function,
	    -analysis_event_id 	     => $event_id,
	});
    }

    $return = $genome_in;
    $return = $return->prepare_for_return();

    #END call_features_insertion_sequences
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_insertion_sequences:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_rRNA_SEED

  $genome_out = $obj->call_features_rRNA_SEED($genome_in, $types)

=over 4




=item Description

Given a genome typed object, find instances of ribosomal RNAs in
the genome.
The types parameter is used to select the types of RNAs to
call. It is a list of strings where each value is one of
   "5S"
   "SSU"
   "LSU"
or "ALL" to choose all available rRNA types.
=back

=cut

sub call_features_rRNA_SEED
{
    my $self = shift;
    my($genome_in, $types) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($types) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"types\" (value was \"$types\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_rRNA_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_rRNA_SEED
    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my %types = map { lc($_) => 1 } @$types;

    my @type_args;

    if ($types{all} || @$types == 0)
    {
	# Don't need to set type arg since tool defaults to calling all.
    }
    else
    {
	for my $type (qw(5S SSU LSU))
	{
	    if ($types{lc($type)})
	    {
		push(@type_args, "-$type");
	    }
	}
    }

    my @cmd = ("rast_call_rRNAs", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id}, @type_args);
    $ctx->stderr->log_cmd(@cmd);
    my $ok = run(\@cmd, $ctx->stderr->redirect);
    if (!$ok)
    {
	die "error calling rRNAs: $?\non command @cmd\n" . $ctx->stderr->text_value;
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    #END call_features_rRNA_SEED
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_rRNA_SEED:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_tRNA_trnascan

  $genome_out = $obj->call_features_tRNA_trnascan($genome_in)

=over 4




=item Description

Given a genome typed object, find instances of tRNAs in
the genome.
=back

=cut

sub call_features_tRNA_trnascan
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_tRNA_trnascan:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_tRNA_trnascan
    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_tRNAs_trnascan", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id});

    $ctx->stderr->log_cmd(@cmd);

    my $ok = run(\@cmd, $ctx->stderr->redirect);

    if (!$ok)
    {
	die "error calling tRNAs: $?\non command @cmd\n" . $ctx->stderr->text_value;
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    #END call_features_tRNA_trnascan
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_tRNA_trnascan:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_RNAs

  $genome_out = $obj->call_RNAs($genome_in)

=over 4




=item Description

Given a genome typed object, find instances of all RNAs we currently
have support for detecting.
=back

=cut

sub call_RNAs
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_RNAs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_RNAs

    $genome_out = $self->call_features_rRNA_SEED($genome_in, 'ALL');
    $genome_out = $self->call_features_tRNA_trnascan($genome_out);

    #END call_RNAs
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_RNAs:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_vigor4

  $return = $obj->call_features_vigor4($genomeTO, $params)

=over 4




=item Description


=back

=cut

sub call_features_vigor4
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_vigor4:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_vigor4

    #
    # Invoke vigor4 to annotate viral genome.
    #
    # Write contigs as fasta file.
    # Invoke vigor4, using the reference that was passed in as the reference_name parameter;
    # We parse the .pep file that is generated. It contains two feature types.
    # This is a CDS:
    # >NC_045512.1 location=266..13468,13471..21555 codon_start=1 gene="orf1ab" ref_db="covid19" ref_id="YP_009724389.1"
    # This is a mature peptide:
    # >NC_045512.1.1 mat_peptide location=266..805 gene="orf1ab" product="leader protein" ref_db="covid19_orf1ab_mp" ref_id="YP_009725297.1"
    #
    # We create features of type CDS and mat_peptide.
    #

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();

    my $ref_dir = $self->{vigor_reference_db_directory};
    $ref_dir ne '' or die "Vigor reference directory not configured";
    -d $ref_dir or die "Vigor reference directory '$ref_dir' not found on $self->{hostname}";

    my @vigor_params = ("-i", $sequences_file,
			"--reference-database-path", $ref_dir,
			"-d", $params->{reference_name},
			"-o", "vigor_out");

    print STDERR Dumper(\@vigor_params);
    # system("cp", "/home/olson/P3/ncov.pep", "vigor_out.pep");
    my $ok = run(["vigor4", @vigor_params],
		 ">", "vigor4.stdout.txt",
		"2>", "vigor4.stderr.txt");
    if (!$ok)
    {
	print STDERR "Vigor run failed with rc=$?. Stdout:\n";
	if (open(V, "<", "vigor4.stdout.txt"))
	{
	    while (<V>)
	    {
		print $_;
	    }
	}
	close(V);
	print STDERR "Stderr:\n";
	if (open(V, "<", "vigor4.stderr.txt"))
	{
	    while (<V>)
	    {
		print $_;
	    }
	}
	close(V);
    }
    
    my $event = {
	tool_name => "vigor4",
	execution_time => scalar gettimeofday,
	parameters => \@vigor_params,
	hostname => $self->{hostname},
    };
    
    my $idc = IDclient->new($genome_in);
    
    my $event_id = $genome_in->add_analysis_event($event);

    #
    # Parse the generated peptide file. We collect the CDS and mature_peptides, then
    # add features so that we can register the counts.
    #
    if (open(my $pep_fh, "<", "vigor_out.pep"))
    {
	my %features;
	while (my($id, $def, $seq) = read_next_fasta_seq($pep_fh))
	{
	    my $fq = { truncated_begin => 0, truncated_end => 0 };
	    
	    my $type;
	    my $ctg;
	    if ($def =~ s/^mat_peptide\s+//)
	    {
		($ctg) = $id =~ /^(.*)\.[^.]+\.[^.]$/;
		$type = 'mat_peptide';
	    }
	    else
	    {
		($ctg) = $id =~ /^(.*)\.[^.]+$/;
		$type = 'CDS';
	    }
	    
	    if (!$ctg)
	    {
		print STDERR "Falling back to prefix of id for contig name from $id\n";
		($ctg) = $id =~ /^(.*?)\./;
	    }
	    
	    my $feature = {
		quality => $fq,
		type => $type,
		contig => $ctg,
		aa_sequence => $seq,
	    };
	    push(@{$features{$type}}, $feature);
	    
	    while ($def =~ /([^=]+)=((\"([^\"]+)\")|([^\"\s]+))\s*/mg) 
	    {
		my $key = $1;
		my $val = $4 ? $4 : $5;
		
		my @loc;
		
		if ($key eq 'location')
		{
		    $feature->{genbank_feature} = { genbank_type => $type, genbank_location  => $val, values => {}};
		    
		    # location=266..13468,13471..21555
		    for my $ent (split(/,/, $val))
		    {
			if (my($s_frag, $s, $e_frag, $e) = $ent =~ /^(<?)(\d+)\.\.(>?)(\d+)$/)
			{
			    $fq->{truncated_begin} = 1 if $s_frag;
			    $fq->{truncated_end} = 1 if $e_frag;
			    
			    my $len = abs($s - $e) + 1;
			    my $strand = $s < $e ? '+' : '-';
			    push(@loc, [$ctg, $s, $strand, $len]);
			}
			else
			{
			    die "error parsing location '$ent'\n";
			}
		    }
		    $feature->{location} = \@loc;
		}
		else
		{
		    $feature->{$key} = $val;
		}
	    }
            $feature->{product} //= $feature->{gene};
	}
	#print Dumper(\%features);
	
	for my $type (keys %features)
	{
	    my $feats = $features{$type};
	    my $n = @$feats;
	    my $id_type = $type;
	    
	    my $id_prefix = $genome_in->{id};
	    if ($id_prefix =~ /^\d+\.\d+$/)
	    {
		$id_prefix = "fig|$id_prefix";
	    }
	    my $typed_prefix = join(".", $id_prefix, $id_type);
	    
	    my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $n);
	    
	    for my $feature (@$feats)
	    {
		my $id = join(".", $typed_prefix, $cur_id_suffix);
		$cur_id_suffix++;
		
		my $p = {
		    -id		     => $id,
		    -type 	     => $type,
		    -location 	     => $feature->{location},
		    -analysis_event_id 	     => $event_id,
		    -annotator => 'vigor4',
		    -protein_translation => $feature->{aa_sequence},
		    -alias_pairs => [[gene => $feature->{gene}]],
		    -function => $feature->{product},
		    -quality_measure => $feature->{quality},
		    -genbank_feature => $feature->{genbank_feature},
		};
		#die Dumper($p);
		
		$genome_in->add_feature($p);
	    }
	}

    }
    else
    {
	warn "Could not read vigor_out.pep\n";
    }
    $return = $genome_in;
    $return = $return->prepare_for_return();

    #END call_features_vigor4
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_vigor4:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_CDS_glimmer3

  $return = $obj->call_features_CDS_glimmer3($genomeTO, $params)

=over 4




=item Description


=back

=cut

sub call_features_CDS_glimmer3
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_glimmer3:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_glimmer3

    my $genome_in = GenomeTypeObject->initialize($genomeTO);

    #
    # Glimmer does not handle code 25.
    #

    if ($genome_in->{genetic_code} != 11 && $genome_in->{genetic_code} != 4)
    {
	print STDERR "Skipping glimmer for genetic code that is neither 11 nor 4\n";
    }
    else
    {
	my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
	
	my $gparams = { %$params };
	$gparams->{genetic_code} = $genome_in->{genetic_code};
	my $calls;
	
	my $stderr = capture_stderr {
	    $calls = Bio::KBase::GenomeAnnotation::Glimmer::call_genes_with_glimmer($sequences_file, $gparams, $ctx);
	};
	if (!$ctx->stderr->log($stderr))
	{
	    print STDERR $stderr;
	}
	
	unlink($sequences_file);
	
	my $trans_table = SeedUtils::genetic_code($genome_in->{genetic_code});
	
	my $event = {
	    tool_name => "glimmer3",
	    execution_time => scalar gettimeofday,
	    parameters => [ map { join("=", $_, $gparams->{$_}) } sort keys %$gparams ],
	    hostname => $self->{hostname},
	};
	
	my $idc = IDclient->new($genome_in);
	
	my $event_id = $genome_in->add_analysis_event($event);
	my $type = 'CDS';
	my $id_type = $self->{cds_id_type};
	
	my $id_prefix = $genome_in->{id};
	if ($id_prefix =~ /^\d+\.\d+$/)
	{
	    $id_prefix = "fig|$id_prefix";
	}
	my $typed_prefix = join(".", $id_prefix, $id_type);
	
	my $count = @$calls;
	my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);
	
	for my $call (@$calls)
	{
	    my($fid, $contig, $begin, $end, $dna) = @$call;
	    
	    my($strand, $len);
	    my $fix_start = 1;
	    if ($begin < $end)
	    {
		$fix_start = 0 if $begin <= 3;
		
		$strand = '+';
		$len = $end - $begin + 1;
	    }
	    else
	    {
		my $cobj = $genome_in->find_contig($contig);
		if (ref $cobj)
		{
		    my $clen = length($cobj->{dna});
		    $fix_start = 0 if $begin > ($clen - 3);
		}
		
		$strand = '-';
		$len = $begin - $end + 1;
	    }
	    
	    my $loc = [[$contig, $begin, $strand, $len]];
	    
	    my $trans = SeedUtils::translate($dna, $trans_table, $fix_start);
	    $trans =~ s/\*$//;
	    
	    my $id = join(".", $typed_prefix, $cur_id_suffix);
	    $cur_id_suffix++;
	    
	    $genome_in->add_feature({
		-id		     => $id,
		-type 	     => $type,
		-location 	     => $loc,
		-analysis_event_id 	     => $event_id,
		-annotator => 'glimmer3',
		-protein_translation => $trans,
	    });
	}
    }
				
    $return = $genome_in;
    $return = $return->prepare_for_return();

    
    #END call_features_CDS_glimmer3
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_glimmer3:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_CDS_prodigal

  $return = $obj->call_features_CDS_prodigal($genomeTO)

=over 4




=item Description


=back

=cut

sub call_features_CDS_prodigal
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_prodigal:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_prodigal

    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genomeTO));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my $stderr = $ctx->stderr;

    my $id_prefix = $genomeTO->{id};
    print "1 '$id_prefix'\n";
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
       print "2 '$id_prefix'\n";
    }
    print "3 '$id_prefix'\n";

    my @cmd = ("rast_call_CDSs_using_prodigal",
	       "--input", $tmp_in,
	       "--output", $tmp_out,
	       "--id-prefix", $id_prefix,
	       "--id-type", $self->{cds_id_type});
    $stderr->log(join(" ", @cmd));

    my $ok = run(\@cmd, $stderr->redirect);
    
    if (!$ok)
    {
	die "error calling CDSs: $?\non command @cmd\n" . $stderr->text_value . "\n";
    }

    $return = $coder->decode(scalar read_file("" . $tmp_out));
    
    #END call_features_CDS_prodigal
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_prodigal:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_CDS_genemark

  $return = $obj->call_features_CDS_genemark($genomeTO)

=over 4




=item Description


=back

=cut

sub call_features_CDS_genemark
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_genemark:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_genemark

    #
    # Genemark requires a license file installed in the "home directory".
    # We expect a parameter to be passed that defines this directory.
    # We also assume that this home directory is the genemark code release. This
    # simplifies things since we can't freely distribute the code in a KB runtime anyway.
    #
    my $gmark_home = $self->{genemark_home};

    #
    # Set up the environment for our subprocess.
    #
    my %env = (HOME => $gmark_home);

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    my($gc) = $genome_in->compute_contigs_gc();

    $gc = int($gc + 0.5);
    print STDERR "gc=$gc\n";

    #
    # Genemark 3.2.5 supports gc 30-70 inclusive.
    #
    $gc = 30 if $gc < 30;
    $gc = 70 if $gc > 70;
    
    my $code = $genome_in->{genetic_code} || 11;
    if ($code != 11 && $code != 4)
    {
	die "Genetic code must be either 4 or 11 for genemark";
    }

    my $h_file = "heu_${code}_${gc}.mod";
    my $h_path = "$gmark_home/heuristic_mod/$h_file";

    if (! -f $h_path)
    {
	die "Heuristic file $h_path not found";
    }

    my $tmp_in = $genome_in->extract_contig_sequences_to_temp_file();
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("$gmark_home/gmhmmp",
	       "-f" => "G",
	       "-r",
	       "-a",
	       "-g", $code,
	       "-m", $h_path,
	       "-o", "$tmp_out",
	       $tmp_in);

    $ctx->stderr->log(join(" ", @cmd)) if $ctx;

    my @cmds = (\@cmd,
		init => sub {
		    $ENV{$_} = $env{$_} foreach keys %env;
		});

    my $event = {
	tool_name => "genemark",
	execution_time => scalar gettimeofday,
	parameters => \@cmd,
	hostname => $self->{hostname},
    };
    my $event_id = $genome_in->add_analysis_event($event);
    my $type = 'CDS';

    my $ok = run(@cmds, ($ctx ? $ctx->stderr->redirect : ()));

    if (!$ok)
    {
	die "Error running pipeline: @cmd\n" . $ctx->stderr->text_value . "\n" if $ctx;
    }

    my $fh;
    open($fh, "<", $tmp_out) or die "Cannot open $tmp_out: $!";

    #
    # Parse generated GFF.
    #
    my $l;
    while (defined($l = <$fh>))
    {
	last unless $l =~ /^#/;
    }

    my %by_gene;
    while (defined($l))
    {
	next if $l =~ /^\s*$/;
	next if $l =~ /^\#/;

	chomp $l;

	my(@fields) = split(/\t/, $l);

	@fields == 9 or die "Invalid GFF at line $.";
	
	my($ctg, $who, $type, $start, $end, $score, $strand, $frame, $attr) = @fields;

	my @attrs = split(/,\s*/, $attr);
	my %attrs = map { my @a = split(/=/, $_, 2); $a[0] => $a[1] } @attrs;

	my $loc;
	my $len = $end - $start + 1;
	if ($strand eq '+') {
	    $loc = [[ $ctg, $start, '+', $len]];
	} else {
	    $loc = [[ $ctg, $end, '-', $len]];
	}

	my $info = [-type => 'CDS',
		    -id_type => $self->{cds_id_type},
		    -location => $loc,
		    -annotator => $who,
		    -annotation => "Add feature called by $who using $h_file",
		    -analysis_event_id => $event_id,
		    -quality_measure => { genemark_score => $score },
		    ];
	$by_gene{$attrs{gene_id}} = $info;
    } continue {
	$l = <$fh>;
    }

    # print STDERR Dumper(\%by_gene);
   

    close($fh);
    open($fh, "<", $tmp_out) or die "Cannot open $tmp_out: $!";
    my $l;
    while (defined($l = <$fh>))
    {
	last unless $l =~ /^#/;
    }

    #
    # We've read features. Set up for creating the IDs etc, then read the proteins and
    # add features to the GTO.
    #
    
    my $id_prefix = $genome_in->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }
    my $typed_prefix = join(".", $id_prefix, $self->{cds_id_type});

    my $count = int(%by_gene);
    my $idc = IDclient->new($genome_in);
    my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

    my $cur;
    my $cur_type;
    my $cur_prot;
    while (defined($l)) {
	next if $l !~ /^\#\#/;
	    
	if ($l =~ /^\#\#(\S+)\s+(\S+)/)
	{
	    $cur_type = $1;
	    $cur = $2;
	    $cur_prot = '';
	    # print STDERR "Set cur_type=$cur_type  cur=$cur\n";
	}
	elsif ($l =~ /^\#\#end-(\S+)/)
	{
	    # print "End: 1=$1 cur=$cur cur_type=$cur_type\n";
	    if ($1 eq 'Protein' && $cur_type eq 'Protein')
	    {
		my $info = $by_gene{$cur};
		if (!$info)
		{
		    warn "No info for gene $info";
		    next;
		}

		my $id = join(".", $typed_prefix, $cur_id_suffix);
		$cur_id_suffix++;

		$genome_in->add_feature({
		    -id => $id,
		    -protein_translation => $cur_prot,
		    @$info,
		});
	    }
	    undef $cur_type;
	    undef $cur;
	    undef $cur_prot;
	}
	elsif ($l =~ /^\#\#(\S+)/)
	{
	    if ($cur_type eq 'Protein')
	    {
		$cur_prot .= $1;
	    }
	    else
	    {
		print STDERR "Skip '$cur_type' $l\n";
	    }
	}
    } continue {
	$l = <$fh>;
    }
    
    $return = $genome_in;
    $return->prepare_for_return();

    #END call_features_CDS_genemark
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_genemark:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_CDS_phanotate

  $return = $obj->call_features_CDS_phanotate($genomeTO)

=over 4




=item Description


=back

=cut

sub call_features_CDS_phanotate
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_phanotate:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_phanotate

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();

    my $trans_table = SeedUtils::genetic_code($genome_in->{genetic_code});

    my $phanotate = which("phanotate");
    if (!$phanotate)
    {
	die "Could not find phanotate in path\n";
    }

    my $version = "Unknown";

    #
    # New phanotate is in pip so we need to use pip to get the version.
    #
    if (open(my $fh, "-|", "pip3", "show", "phanotate"))
    {
	while (<$fh>)
	{
	    if (/^Version:\s+(\S+)/)
	    {
		$version = $1;
	    }
	}
	close($fh);
    }

    my $event = {
	tool_name => "$phanotate version $version",
	execution_time => scalar gettimeofday,
	parameters => $sequences_file,
	hostname => $self->{hostname},
    };

    my $idc = IDclient->new($genome_in);
    
    my $event_id = $genome_in->add_analysis_event($event);
    my $type = 'CDS';
    my $id_type = $self->{cds_id_type};

    my $id_prefix = $genome_in->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }

    my @calls;
    my $typed_prefix = join(".", $id_prefix, $id_type);
    if (open(my $fh, "-|", $phanotate, $sequences_file))
    {
	while (<$fh>)
	{
	    chomp;
	    next if /^#/;
	    my($start, $stop, $strand, $contig, $score) = split(/\t/);

	    push(@calls, [$start, $stop, $strand, $contig, $score]);
	}
	close($fh);
    }
    else
    {
	die "Error starting phanotate $sequences_file: $!";
    }

    #
    # Determine contig lengths for fixing fuzzy coordinates.
    #
    my %contig_length;
    for my $contig ($genome_in->contigs)
    {
	$contig_length{$contig->{id}} = length($contig->{dna});
    }

    my $count = @calls;
    my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

    for my $call (@calls)
    {
	my($start, $stop, $strand, $contig, $score) = @$call;

	my $fix_start = 1;

	$start = _fix_fuzzy_coordinate($start, $contig_length{$contig});
	$stop = _fix_fuzzy_coordinate($stop, $contig_length{$contig});

	my($left, $right) = $strand eq '+' ? ($start, $stop) : ($stop, $start);

	my $len = $right - $left + 1;
	
	my $cobj = $genome_in->find_contig($contig);
	if (ref $cobj)
	{
	    my $loc = [[$contig, $start, $strand, $len]];
	    
	    my $dna = $genome_in->get_feature_dna({ location => $loc });
	    my $trans = SeedUtils::translate($dna, $trans_table, $fix_start);
	    $trans =~ s/\*$//;
	    
	    my $id = join(".", $typed_prefix, $cur_id_suffix);
	    $cur_id_suffix++;
	    
	    $genome_in->add_feature({
		-id		     => $id,
		-type 	     => $type,
		-location 	     => $loc,
		-analysis_event_id 	     => $event_id,
		-annotator => 'phanotate',
		-protein_translation => $trans,
	    });
	}
    }

    unlink($sequences_file);

    $return = $genome_in;
    $return = $return->prepare_for_return();


    #END call_features_CDS_phanotate
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_phanotate:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 prune_invalid_CDS_features

  $genome_out = $obj->prune_invalid_CDS_features($genome_in, $params)

=over 4




=item Description


=back

=cut

sub prune_invalid_CDS_features
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to prune_invalid_CDS_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN prune_invalid_CDS_features

    #
    # Find and remove any feature whose DNA is mostly a single character.
    # We do this only on small features since this mostly a problem
    # with assembler-generated tiny contigs.
    #

    $genome_in = GenomeTypeObject->initialize($genome_in);

    #
    # Compute contig lengths to prune based on that.
    #
    my $min_contig_length = $params->{minimum_contig_length} // 0;
    my $max_homopolymer_frequency = $params->{max_homopolymer_frequency} // 0.9;
    my %contig_length;
    for my $contig ($genome_in->contigs)
    {
	$contig_length{$contig->{id}} = length($contig->{dna});
    }

    my %examine = (CDS => 1, peg => 1);
    my @delete;
    for my $f (grep { $examine{$_->{type}} } @{$genome_in->features})
    {
	my($contig, $beg, $strand, $len) = @{$f->{location}->[0]};
	next if $len > 300;
	if ($contig_length{$contig} < $min_contig_length)
	{
	    printf STDERR "Rejecting $f->{id} due contig length $contig_length{$contig} < $min_contig_length\n";
	    push(@delete, $f->{id});
	    next;
	}
	
	my $dna = $genome_in->get_feature_dna($f);
	#
	# compute freqs of each char.
	#
	my $fa = ($dna =~ tr/aA//) / $len;
	my $fc = ($dna =~ tr/cC//) / $len;
	my $fg = ($dna =~ tr/gG//) / $len;
	my $ft = ($dna =~ tr/tT//) / $len;
	if ($fa > $max_homopolymer_frequency ||
	    $fc > $max_homopolymer_frequency ||
	    $fg > $max_homopolymer_frequency ||
	    $ft > $max_homopolymer_frequency)
	{
	    printf STDERR "Rejecting $f->{id} due to freqs A=%.2f C=%2.f G=%.2f T=%.2f\n",
	    $fa, $fc, $fg, $ft;
	    push(@delete, $f->{id});
	}
    }

    for my $fid (@delete)
    {
	$genome_in->delete_feature($fid);
    }
    
    $genome_out = $genome_in->prepare_for_return();

    #END prune_invalid_CDS_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to prune_invalid_CDS_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_CDS_SEED_projection

  $return = $obj->call_features_CDS_SEED_projection($genomeTO, $params)

=over 4




=item Description


=back

=cut

sub call_features_CDS_SEED_projection
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_SEED_projection:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_SEED_projection

    #
    # Do some validity checking.
    #

    #
    # Default to coreseed.
    #
    $params->{reference_database} //= 'core';
    my $ref_db = $params->{reference_database};

    if ($ref_db eq 'PATRIC')
    {
	die "PATRIC reference database not yet supported\n";
    }

    if ($params->{reference_id} eq '')
    {
	die "No reference ID provided\n";
    }

    if ($ref_db !~ /^http/)
    {
	my $where = SeedURLs::url($ref_db);
	if (!$where)
	{
	    die "Could not find reference database '$ref_db'\n";
	}
    }

    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genomeTO));

    my @cmd = ("fast_project",
	       "--input", $tmp_in,
	       "-s", $ref_db,
	       "-r", $params->{reference_id},
	       (defined($params->{kmer_size}) ? ("--kmersize", $params->{kmer_size}) : ()),
	       "--format", "gto");

    $ctx->stderr->log_cmd(@cmd);
    my $genomeOut_json;
    my $ok = run(\@cmd,
		 '>', \$genomeOut_json,
		 $ctx->stderr->redirect);

    if (!$ok)
    {
	die "error calling fast_project: $?\non command @cmd";
    }

    $return = $coder->decode($genomeOut_json);


    #END call_features_CDS_SEED_projection
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_SEED_projection:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_CDS_FragGeneScan

  $return = $obj->call_features_CDS_FragGeneScan($genomeTO)

=over 4




=item Description


=back

=cut

sub call_features_CDS_FragGeneScan
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_FragGeneScan:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_FragGeneScan
    die "Not implemented";
    #END call_features_CDS_FragGeneScan
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_FragGeneScan:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_repeat_region_SEED

  $genome_out = $obj->call_features_repeat_region_SEED($genome_in, $params)

=over 4




=item Description


=back

=cut

sub call_features_repeat_region_SEED
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_repeat_region_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_repeat_region_SEED

    $genome_in = GenomeTypeObject->initialize($genome_in);

    #
    # We skip this step if we have too many contigs.
    #
    if ($genome_in->n_contigs > 10)
    {
	print STDERR "Skipping call_features_repeat_region_SEED because contig count " . $genome_in->n_contigs . " is too large\n";
    }
    else
    {
	my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
	my $output_file = File::Temp->new();
	
	my @opts;
	push(@opts, "-i", $params->{min_identity}) if defined($params->{min_identity});
	push(@opts, "-l", $params->{min_length}) if defined($params->{min_length});
	
	my $stderr = $ctx->stderr;
	
	my @cmd = ("svr_big_repeats", @opts);
	
	$stderr->log_cmd(@cmd);
	
	my $tmpdir = File::Temp->newdir(undef, CLEANUP => 1);
	
	print STDERR "seq file $sequences_file\n";
	my $ok = run(\@cmd,
		     init => sub { print STDERR "init chdir to $tmpdir\n"; chdir($tmpdir) or die $!; },
		     "<", $sequences_file,
		     "|",
		     ["svr_condense_repeats"],
		     ">", $output_file,
		     $stderr->redirect,
		    );
	
	unlink($sequences_file);
	
	if (!$ok)
	{
	    die "Error running svr_big_repeats: @cmd\n" . $stderr->text_value . "\n";
	}
	
	close($output_file);
	my($res_fh);
	open($res_fh, "<", $output_file) or die "Cannot open svr_big_repeats file $output_file: $!";
	
	my $event = {
	    tool_name => "svr_big_repeats | svr_condense_repeats",
	    execution_time => scalar gettimeofday,
	    parameters => \@opts,
	    hostname => $self->{hostname},
	};
	
	my $event_id = $genome_in->add_analysis_event($event);
	
	# olson@bio-data-1:~/FIGdisk/dist/releases/dev2$ svr_condense_repeats < r
	#    NC_000913       15377   16741   99.71
	#    NC_000913       19796   20564   98.83
	
	my $type = 'repeat';
	my $function = 'repeat region';
	
	while(<$res_fh>)
	{
	    chomp;
	    my($contig, $left, $right, $iden) = split(/\t/);
	    
	    next unless $left =~ /^\d+$/ && $right =~ /^\d+$/;
	    
	    my $confidence = $iden / 100.0;
	    
	    my $quality = {
		existence_confidence => $confidence,
	    };
	    
	    my $len = 1 + $right - $left;
	    my $loc = [[$contig, $left, '+', $len]];
	    $genome_in->add_feature({
		-id_prefix 	     => $genome_in->{id},
		-type 	     => $type,
		-location 	     => $loc,
		-function 	     => $function,
		-analysis_event_id 	     => $event_id,
		-quality_measure => $quality,
	    });
	}
    }
    
    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();

    #END call_features_repeat_region_SEED
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_repeat_region_SEED:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_prophage_phispy

  $genome_out = $obj->call_features_prophage_phispy($genome_in)

=over 4




=item Description


=back

=cut

sub call_features_prophage_phispy
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_prophage_phispy:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_prophage_phispy
    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_prophage_using_phispy", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id});

    my $stderr = $ctx->stderr;
    $stderr->log_cmd(@cmd);

    my $ok = run(\@cmd, $stderr->redirect);

    if (!$ok)
    {
	die "error calling prophages: $?\non command @cmd\n" . $stderr->text_value . "\n";
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    #END call_features_prophage_phispy
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_prophage_phispy:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_scan_for_matches

  $genome_out = $obj->call_features_scan_for_matches($genome_in, $pattern, $feature_type)

=over 4




=item Description


=back

=cut

sub call_features_scan_for_matches
{
    my $self = shift;
    my($genome_in, $pattern, $feature_type) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (!ref($pattern)) or push(@_bad_arguments, "Invalid type for argument \"pattern\" (value was \"$pattern\")");
    (!ref($feature_type)) or push(@_bad_arguments, "Invalid type for argument \"feature_type\" (value was \"$feature_type\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_scan_for_matches:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_scan_for_matches
    die "Not implemented";
    #END call_features_scan_for_matches
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_scan_for_matches:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_assembly_gap

  $genome_out = $obj->call_features_assembly_gap($genome_in, $params)

=over 4




=item Description

Given a genome typed object, call gap features.                
Gaps are known regions in the contig where the nucleotide sequence is not known                
but where there is evidence that a run of DNA does exist joining the sequenced                        
data on either side of the gap.    
        
Gaps are currently called using one of two methods. Genomes that originated as                
genbank files may have a CONTIGS entry that defines the contig and gap regions.                        
Genomes that do not have a CONTIGS entry are scanned for runs of "n" characters.
=back

=cut

sub call_features_assembly_gap
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_assembly_gap:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_assembly_gap

    #
    # If this genome already has assembly_gap features, we will just return the genome
    # unchanged.
    #
    # Otherwise, we will scan the contigs looking for runs.
    #

    for my $f (@{$genome_in->{features}})
    {
	if ($f->{type} eq 'assembly_gap')
	{
	    print STDERR "Genome already has assembly gap features\n";
	    $genome_out = $genome_in;
	    goto done;
	}
    }
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $event = {
	tool_name => "call_features_assembly_gap",
	execution_time => scalar gettimeofday,
	parameters => [ map { join("=", $_, $params->{$_}) } sort keys %$params ],
	hostname => $self->{hostname},
    };
    my $idc = IDclient->new($genome_in);
    
    my $event_id = $genome_in->add_analysis_event($event);
    my $type = 'assembly_gap';
    my $id_prefix = $genome_in->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }
    my $typed_prefix = join(".", $id_prefix, $type);

    my $split = $params->{min_gap_length} // 4;
    my $longrepeat = $params->{monopolymer_repeat_length} // 50;

    for my $ctg (@{$genome_in->{contigs}})
    {
	my $dna = $ctg->{dna};
	
	while ($dna =~ m/([nbdhvrykmswx]{$split,}|a{$longrepeat,}|c{$longrepeat,}|g{$longrepeat,}|t{$longrepeat,})/gio)
	{
	    my $len = length($1);
	    my $end = pos($dna);
	    my $start = $end - $len + 1;

	    my $loc = [[$ctg->{id}, $start, '+', $len]];
	    my $function = "assembly gap";

	    $genome_in->add_feature({
		-id_prefix 	   => $genome_in->{id},
		-type		   => $type,
		-location 	   => $loc,
		-function 	   => $function,
		-analysis_event_id => $event_id,
	    });
	}
    }
    
    $genome_out= $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
 done:
    #END call_features_assembly_gap
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_assembly_gap:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 split_gap_spanning_features

  $genome_out = $obj->split_gap_spanning_features($genome_in, $split_gap_spanning_features_params)

=over 4




=item Description


=back

=cut

sub split_gap_spanning_features
{
    my $self = shift;
    my($genome_in, $split_gap_spanning_features_params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($split_gap_spanning_features_params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"split_gap_spanning_features_params\" (value was \"$split_gap_spanning_features_params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to split_gap_spanning_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN split_gap_spanning_features

    #
    # Collect the assembly gap features in this genome, sorted in order by start.
    # For each other feature, find overlaps with the gaps. Split them the range of the gap.
    # If we split proteins, remove the translation. We will re-translate in a later stage.
    #

    $genome_in = GenomeTypeObject->initialize($genome_in);

    my %gaps_by_contig;
    for my $gap (grep { $_->{type} eq 'assembly_gap' } @{$genome_in->{features}})
    {
	my $loc = $gap->{location}->[0];
	my($ctg, $start, $strand, $len) = @$loc;
	my $end = $start + $len - 1;
	push @{$gaps_by_contig{$ctg}}, [$gap->{id}, $start, $end, $len];
    }

    my %type_to_check = (CDS => 1, peg => 1, rna => 1 );

    my $id_server = IDclient->new($genome_in);

    for my $feature (grep { $type_to_check{$_->{type}} } @{$genome_in->{features}})
    {
	my $loc = $feature->{location};
	#
	# Don't try to split multi-location features at this point.
	#
	next if @$loc > 1;
	$loc = $loc->[0];
	my($ctg, $fstart, $strand, $flen) = @$loc;

	my($fleft, $fright);
	if ($strand eq '+')
	{
	    $fleft = $fstart;
	    $fright = $fleft + $flen - 1;
	}
	else
	{
	    $fright = $fstart;
	    $fleft = $fstart - $flen + 1;
	}
	    
	my $gaplist = $gaps_by_contig{$loc->[0]};

	for my $gap (@$gaplist)
	{
	    my($gap_id, $gstart, $gend, $glen) = @$gap;
	    
	    if ($gend >= $fleft && $gstart <= $fright)
	    {
		print STDERR "$feature->{id} overlaps gap $gap_id\n";
		print Dumper($feature->{id}, $loc, $gap);

		# Determine gap is an insert or just truncates the feature.

		if ($fright > $gend && $fleft < $gstart)
		{
		    #
		    # Full overlap. Create two features.
		    #
		    my $f1left = $fleft;
		    my $f1right = $gstart - 1;
		    my $f2left = $gend + 1;
		    my $f2right = $fright;

		    my $f1len  = $f1right - $f1left + 1;
		    my $f2len  = $f2right - $f2left + 1;

		    #
		    # Keep the downstream truncation in-frame
		    #
		    
		    if ($strand eq '+')
		    {
			my $offset = $f2len % 3;
			$f2left += $offset;
			$f2len  = $f2right - $f2left + 1;
		    }
		    else
		    {
			my $offset = $f1len % 3;
			$f1right -= $offset;
			$f1len  = $f2right - $f2left + 1;
		    }
		    

		    my $fleft = Clone::clone($feature);
		    my $id = _allocate_new_feature_id($id_server, $feature->{id}, $fleft->{type});
		    $fleft->{id} = $id;
		    $fleft->{location} = [ _ends_to_location($ctg, $strand, $f1left, $f1right) ];
		    delete $fleft->{protein_translation};
		    my $where_left = ($strand eq '+') ? "truncated_end" : "truncated_begin";
		    $fleft->{$where_left} = 1;
		    
		    my $fright = Clone::clone($feature);
		    my $id = _allocate_new_feature_id($id_server, $feature->{id}, $fright->{type});
		    $fright->{id} = $id;
		    $fright->{location} = [ _ends_to_location($ctg, $strand, $f2left, $f2right) ];
		    delete $fright->{protein_translation};
		    my $where_right = ($strand eq '+') ? "truncated_begin" : "truncated_end";
		    $fright->{$where_right} = 1;

		    push(@{$genome_in->{features}}, $fleft, $fright);
		    print Dumper($fleft, $fright);
		}
		else
		{
		    print "   partial overlap\n";
		}
	    }
	}
    }
    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    #END split_gap_spanning_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to split_gap_spanning_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 translate_untranslated_proteins

  $genome_out = $obj->translate_untranslated_proteins($genome_in)

=over 4




=item Description


=back

=cut

sub translate_untranslated_proteins
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to translate_untranslated_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN translate_untranslated_proteins

    #
    # For every protein that doesn't have a translation, attempt to create one.
    #

    $genome_in = GenomeTypeObject->initialize($genome_in);
    my $code = SeedUtils::genetic_code($genome_in->{genetic_code});

    for my $feature (grep { ($_->{type} eq 'CDS' || $_->{type} eq 'peg') && !$_->{protein_translation} } @{$genome_in->{features}})
    {
	my $dna = $genome_in->get_feature_dna($feature->{id});
	# Don't fix start char if we know we have a truncated start.
	my $trans = SeedUtils::translate($dna, $code, !($feature->{truncated_begin}));
	$trans =~ s/\*$//;
	print STDERR "Translate $feature->{id} $dna => $trans\n";
	$feature->{protein_translation} = $trans;
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    #END translate_untranslated_proteins
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to translate_untranslated_proteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 annotate_proteins_similarity

  $return = $obj->annotate_proteins_similarity($genomeTO, $params)

=over 4




=item Description

Annotate based on similarity to annotation databases.
=back

=cut

sub annotate_proteins_similarity
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins_similarity:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_proteins_similarity

    my $coder = _get_coder();
    my $stderr = $ctx->stderr;
    
    my $dir = $self->{nr_annotation_directory};
    if (! -d $dir)
    {
	warn "NR annotation directory '$dir' empty or missing\n";
	$return = $genomeTO;
    }
    else
    {
	my $genomeTO_json = $coder->encode($genomeTO);

	my $genomeOut_json;

	my $tmp = File::Temp->new();
	print $tmp $genomeTO_json;
	close($tmp);

	print STDERR "Starting $dir $tmp\n";
	my @params = ("--nr-dir", $dir,
		      "--input", $tmp);

	if ($params->{annotate_null_only})
	{
	    push(@params, "-N");
	}
	elsif ($params->{annotate_hypothetical_only})
	{
	    push(@params, "-H");
	}
	
	my @cmd = ('rast_annotate_proteins_similarity', @params);
	my $ok = run(\@cmd,
		     '>', \$genomeOut_json,
		     $stderr->redirect);

	undef $tmp;
	undef $genomeTO;

	if ($ok) {
	    $return = $coder->decode($genomeOut_json);
	} else {
	    die "rast_annotate_proteins_similarity failed ($?):\n" . $stderr->text_value;
	}
    }

    #END annotate_proteins_similarity
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins_similarity:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 annotate_proteins_phage

  $return = $obj->annotate_proteins_phage($genomeTO, $params)

=over 4




=item Description

Annotate based on similarity to the phage annotation daatabase.
=back

=cut

sub annotate_proteins_phage
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins_phage:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_proteins_phage

    my $coder = _get_coder();
    my $stderr = $ctx->stderr;
    
    my $files = $self->{phage_annotation_files};
    if (!ref($files) || @$files == 0)
    {
	warn "No phage annotation files found\n";
	$return = $genomeTO;
    }
    else
    {
	my $tmpcur = File::Temp->new();
	print $tmpcur $coder->encode($genomeTO);
	close($tmpcur);

	undef $genomeTO;

	for my $pfile (@$files)
	{
	    my $tmpout = File::Temp->new();
	    close($tmpout);
		
	    print STDERR "Starting phage annotation $pfile < $tmpcur > $tmpout\n";

	    my @params = ("--nr-file", $pfile, 
			  '--input', $tmpcur,
			  '--output', $tmpout);
	    
	    if ($params->{annotate_null_only})
	    {
		push(@params, "-N");
	    }
	    elsif ($params->{annotate_hypothetical_only})
	    {
		push(@params, "-H");
	    }


	    my @cmd = ('rast_annotate_proteins_similarity', @params);
	    
	    my $ok = run(\@cmd,
			 $stderr->redirect);

	    if (!$ok)
	    {
		die "rast_annotate_proteins_similarity failed during phage annotation ($?) cmd=@cmd:\n" . $stderr->text_value;
	    }

	    undef $tmpcur;
	    $tmpcur = $tmpout;
	}

	$return = $coder->decode(scalar read_file("" . $tmpcur));
    }

    #END annotate_proteins_phage
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins_phage:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 annotate_proteins_kmer_v1

  $return = $obj->annotate_proteins_kmer_v1($genomeTO, $params)

=over 4




=item Description


=back

=cut

sub annotate_proteins_kmer_v1
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_proteins_kmer_v1

    my $n_proteins_per_call = 100;

    my $genome_in = GenomeTypeObject->initialize($genomeTO);

    my $kmer_service = Bio::KBase::KmerAnnotationByFigfam::Client->new($self->{kmer_service_url});


    if (!defined($params->{dataset_name}))
    {
	$params->{dataset_name} = $kmer_service->get_default_dataset_name();
    }

    if (!defined($params->{kmer_size}))
    {
	$params->{kmer_size} = 8;
    }

    warn "Annotating using dataset $params->{dataset_name} kmer_size $params->{kmer_size}\n";
    my $event = {
	tool_name => "KmerAnnotationByFigfam",
	execution_time => scalar gettimeofday,
	parameters => [ map { join("=", $_, $params->{$_}) } sort keys %$params ],
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);

    my @proteins;

    my $do_anno = sub {
	my($proteins) = @_;
	my $res = $kmer_service->annotate_proteins($proteins, $params);
	for my $hit (@$res)
	{
	    my($id, $func, $otu, $score, $nonover, $over, $details) = @$hit;
	    if ($func)
	    {
		$genome_in->update_function("annotate_proteins_kmer_v1", $id, $func, $event_id);
		my $feature = $genome_in->find_feature($id);
		if (ref($feature) eq 'HASH')
		{
		    $feature->{ quality }->{ hit_count } = $score;
		}
	    }
	}
    };

    for my $feature ($genome_in->features)
    {
	my $trans = $feature->{protein_translation};
	next unless $trans;

	if ($params->{annotate_null_only})
	{
	    my $f = $feature->{function};
	    if ($f ne '')
	    {
		next;
	    }
	}
	elsif ($params->{annotate_hypothetical_only})
	{
	    my $f = $feature->{function};
	    if ($f ne '' && $f !~ /^\s*hypothetical\s+protein\s*$/i)
	    {
		next;
	    }
	}

	push(@proteins, [$feature->{id}, $trans]);
	if (@proteins >= $n_proteins_per_call)
	{
	    $do_anno->(\@proteins);
	    @proteins = ();
	}
    }
    $do_anno->(\@proteins) if @proteins;
    
    $return = $genome_in->prepare_for_return();

    #END annotate_proteins_kmer_v1
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 annotate_proteins_kmer_v2

  $genome_out = $obj->annotate_proteins_kmer_v2($genome_in, $params)

=over 4




=item Description


=back

=cut

sub annotate_proteins_kmer_v2
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_proteins_kmer_v2

    $genome_in = GenomeTypeObject->initialize($genome_in);

    my $filter;
    if ($params->{annotate_null_only})
    {
	$filter = sub {
	    my($feature) = @_;
	    my $f = $feature->{function};
	    if ($f eq '')
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	};
    }
    elsif ($params->{annotate_hypothetical_only})
    {
	$filter = sub {
	    my($feature) = @_;
	    my $f = $feature->{function};
	    if ($f eq '' || $f =~ /^\s*hypothetical\s+protein\s*$/i)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	};
    }
    my $sequences_file = $genome_in->extract_protein_sequences_to_temp_file($filter);
    my $output_file = File::Temp->new();

    my $min_hits = 5;
    my $max_gap = 200;

    if (defined($params->{min_hits}))
    {
	$min_hits = $params->{min_hits};
    }
	
    if (defined($params->{max_gap}))
    {
	$max_gap = $params->{max_gap};
    }

    my $data_dir = $self->{kmer_v2_data_directory};
    my @service_url;
    
    if ($params->{kmer_data_directory})
    {
	$data_dir = $params->{kmer_data_directory};
	-d $data_dir or die "Kmer data directory $data_dir does not exist";
	if ($params->{kmer_service_url})
	{
	    @service_url = ("-u", $params->{kmer_service_url});
	}
    }
    else
    {
	#
	# If we're configured with a families server and it is configured
	# with our v2 kmers, use it for function assignment.
	#
	if (($self->{patric_annotate_families_kmers} eq $self->{kmer_v2_data_directory}) &&
	    $self->{patric_annotate_families_url})
	{
	    @service_url = ("-u", $self->{patric_annotate_families_url} . "/query");
	}
    }
    my @params = ("-a", "-g", $max_gap, "-m", $min_hits, "-d", $data_dir, @service_url);

    my @cmd = ("kmer_search", @params);
    $ctx->stderr->log_cmd(@cmd);
    my $ok = run(\@cmd,
		 "<", $sequences_file,
		 ">", $output_file,
		 $ctx->stderr->redirect);
    
    unlink($sequences_file);

    if (!$ok)
    {
	die "Error running kmer_search ($?): @cmd\n" . $ctx->stderr->text_value;
    }

    close($output_file);
    my($res_fh);
    open($res_fh, "<", $output_file) or die "Cannot open kmer_search output file $output_file: $!";

    my $event = {
	tool_name => "kmer_search",
	execution_time => scalar gettimeofday,
	parameters => \@params,
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);
    
    while(<$res_fh>)
    {
	chomp;
	#  my($fid, $function) = split(/\t/);
	my($fid, $function, $hits, $hitsW) = split(/\t/);
	$genome_in->update_function("annotate_proteins_kmer_v2", $fid, $function, $event_id);

	my $feature = $genome_in->find_feature($fid);
	if (ref($feature) eq 'HASH')
	{
	    $feature->{ quality }->{ hit_count } = $hits;
	    $feature->{ quality }->{ weighted_hit_count } = $hitsW;
	}
    }

    $genome_out = $genome_in->prepare_for_return();

    #END annotate_proteins_kmer_v2
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 resolve_overlapping_features

  $genome_out = $obj->resolve_overlapping_features($genome_in, $params)

=over 4




=item Description


=back

=cut

sub resolve_overlapping_features
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to resolve_overlapping_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN resolve_overlapping_features

    $genome_in = GenomeTypeObject->initialize($genome_in);

    my $stderr;

    $ctx->stderr->log_cmd("resolve_overlapping_features", Dumper($params));
    $stderr = capture_stderr {
	$genome_out = overlap_resolution::resolve_overlapping_features($genome_in, $params);
    };
    if (!$ctx->stderr->log($stderr))
    {
	print STDERR $stderr;
    }

    $genome_out = $genome_out->prepare_for_return();

    #END resolve_overlapping_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to resolve_overlapping_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 propagate_genbank_feature_metadata

  $genome_out = $obj->propagate_genbank_feature_metadata($genome_in, $params)

=over 4




=item Description


=back

=cut

sub propagate_genbank_feature_metadata
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to propagate_genbank_feature_metadata:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN propagate_genbank_feature_metadata

    $genome_in = GenomeTypeObject->initialize($genome_in);

    my $stderr;

    $ctx->stderr->log_cmd("propagate_genbank_feature_metadata", Dumper($params));
    $stderr = capture_stderr {
	$genome_out = PropagateGBMetadata::propagate_gb_metadata($genome_in, $params, $ctx->user_id);

    };
    if (!$ctx->stderr->log($stderr))
    {
	print STDERR $stderr;
    }

    $genome_out = $genome_out->prepare_for_return();

    #END propagate_genbank_feature_metadata
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to propagate_genbank_feature_metadata:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 call_features_ProtoCDS_kmer_v1

  $return = $obj->call_features_ProtoCDS_kmer_v1($genomeTO, $params)

=over 4




=item Description


=back

=cut

sub call_features_ProtoCDS_kmer_v1
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_ProtoCDS_kmer_v1

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    
    my $min_hits = 5;
    my $max_gap = 200;

    if (defined($params->{min_hits}))
    {
	$min_hits = $params->{min_hits};
    }
	
    if (defined($params->{max_gap}))
    {
	$max_gap = $params->{max_gap};
    }

    my $event = {
	tool_name => "KmerAnnotationByFigfam",
	execution_time => scalar gettimeofday,
	parameters => [ map { join("=", $_, $params->{$_}) } sort keys %$params ],
	hostname => $self->{hostname},
    };

    my $idc = IDclient->new($genome_in);
    
    my $event_id = $genome_in->add_analysis_event($event);

    my $type = 'protoCDS';
    my $id_prefix = $genome_in->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }
    my $typed_prefix = join(".", $id_prefix, $type);

    my $kmer_service = Bio::KBase::KmerAnnotationByFigfam::Client->new($self->{kmer_service_url});
    if (!defined($params->{dataset_name}))
    {
	$params->{dataset_name} = $kmer_service->get_default_dataset_name();
    }

    if (!defined($params->{kmer_size}))
    {
	$params->{kmer_size} = 8;
    }
    
    for my $ctg ($genome_in->contigs)
    {
	my $hits = $kmer_service->call_genes_in_dna([[$ctg->{id}, $ctg->{dna}]], $params);
	# print STDERR Dumper($hits);

	my $count = @$hits;
	my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

	for my $hit (@$hits)
	{
	    my($nhits, $id, $begin, $end, $function, $otu) = @$hit;
	    my $quality = {
		existence_confidence => 0.5,
		hit_count => 0 + $nhits,
	    };
	    my($strand, $len);
	    if ($begin < $end)
	    {
		$strand = '+';
		$len = $end - $begin + 1;
	    }
	    else
	    {
		$strand = '-';
		$len = $begin - $end + 1;
	    }
	    
	    my $fid = join(".", $typed_prefix, $cur_id_suffix);
	    $cur_id_suffix++;

	    my $loc = [[$id, $begin, $strand, $len]];
	    $genome_in->add_feature({
		-id		     => $fid,
		-type 	     => $type,
		-location 	     => $loc,
		-function 	     => $function,
		-analysis_event_id 	     => $event_id,
		-quality_measure => $quality,
	    });
	}
    }

    $return = $genome_in;
    $return = $return->prepare_for_return();

    #END call_features_ProtoCDS_kmer_v1
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_ProtoCDS_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_ProtoCDS_kmer_v2

  $genome_out = $obj->call_features_ProtoCDS_kmer_v2($genome_in, $params)

=over 4




=item Description


=back

=cut

sub call_features_ProtoCDS_kmer_v2
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_ProtoCDS_kmer_v2

    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
    my $output_file = File::Temp->new();

    my $min_hits = 5;
    my $max_gap = 200;

    if (defined($params->{min_hits}))
    {
	$min_hits = $params->{min_hits};
    }
	
    if (defined($params->{max_gap}))
    {
	$max_gap = $params->{max_gap};
    }

    my @params = ("-g", $max_gap, "-m", $min_hits, "-d", $self->{kmer_v2_data_directory});

    my @cmd = ("kmer_search", @params);
    my $ok = run(\@cmd,
		 "<", $sequences_file,
		 ">", $output_file);

    unlink($sequences_file);

    if (!$ok)
    {
	die "Error running kmer_search: @cmd\n";
    }

    close($output_file);
    my($res_fh);
    open($res_fh, "<", $output_file) or die "Cannot open kmer_search output file $output_file: $!";

    my $event = {
	tool_name => "kmer_search",
	execution_time => scalar gettimeofday,
	parameters => \@params,
	hostname => $self->{hostname},
    };

    my $idc = IDclient->new($genome_in);
    my $event_id = $genome_in->add_analysis_event($event);

#	NC_002952       2902278 2902407 -       1       37      LSU ribosomal protein L34p      116.145798
#	Contig
#	Left
#	Right
#	Strand
#	Frame
#	Score-1
#	Function
#	Score-2
#    

    my $type = 'protoCDS';
    my $id_prefix = $genome_in->{id};
    if ($id_prefix =~ /^\d+\.\d+$/)
    {
	$id_prefix = "fig|$id_prefix";
    }
    my $typed_prefix = join(".", $id_prefix, $type);

    my @hits;

    while(<$res_fh>)
    {
	chomp;
	my($contig, $left, $right, $strand, $frame, $hit_count, $function, $weighted_hit_count) = split(/\t/, $_);

	next unless $left =~ /^\d+$/ && $right =~ /^\d+$/;

	push(@hits, $_);
    }
    close($res_fh);

    my $count = @hits;
    my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

    for my $hit (@hits)
    {
	my($contig, $left, $right, $strand, $frame, $hit_count, $function, $weighted_hit_count) = split(/\t/, $hit);

	my $confidence = 1 - 0.5 ** ($weighted_hit_count / 3);

	my $quality = {
	    existence_confidence => $confidence,
	    hit_count => 0 + $hit_count,
	    weighted_hit_count => 0 + $weighted_hit_count,
	};

	my $id = join(".", $typed_prefix, $cur_id_suffix);
	$cur_id_suffix++;

	my $begin = 0 + (($strand eq '+') ? $left : $right);
	my $len = 1 + $right - $left;
	my $loc = [[$contig, $begin, $strand, $len]];
	$genome_in->add_feature({
	    -id              => $id,
	    -type 	     => $type,
	    -location 	     => $loc,
	    -function 	     => $function,
	    -analysis_event_id 	     => $event_id,
	    -quality_measure => $quality,
	});
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();

    #END call_features_ProtoCDS_kmer_v2
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_ProtoCDS_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 enumerate_special_protein_databases

  $database_names = $obj->enumerate_special_protein_databases()

=over 4




=item Description


=back

=cut

sub enumerate_special_protein_databases
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($database_names);
    #BEGIN enumerate_special_protein_databases
    
    $database_names = [];
    
    my $dir = $self->{special_protein_dbdir};
    if ($dir && opendir(my $dh, $dir))
    {
	my @list = map { s/\.faa$//; $_ } grep { -f "$dir/$_" && /\.faa$/ } readdir($dh);
	push(@$database_names, @list);
	closedir($dh);
    }
    
    #END enumerate_special_protein_databases
    my @_bad_returns;
    (ref($database_names) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"database_names\" (value was \"$database_names\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to enumerate_special_protein_databases:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($database_names);
}


=head2 compute_special_proteins

  $results = $obj->compute_special_proteins($genome_in, $database_names)

=over 4




=item Description


=back

=cut

sub compute_special_proteins
{
    my $self = shift;
    my($genome_in, $database_names) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($database_names) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"database_names\" (value was \"$database_names\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to compute_special_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($results);
    #BEGIN compute_special_proteins
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $dir = $self->{special_protein_dbdir};
    
    my $prots = File::Temp->new();
    close($prots);
    $genome_in->write_protein_translations_to_file($prots);
    
    my $out = File::Temp->new();
    
    my @dbargs = map { ("--db", $_) } @$database_names;
    
    if (!$self->{special_protein_cache_db})
    {
	push(@dbargs, "--no-cache");
    }
    
    my @cmd = ('rast_compute_specialty_genes',
	       "--in", $prots,
	       "--db-dir", $dir,
	       @dbargs,
	       $self->{special_protein_cache_db} ? ("--cache-db", $self->{special_protein_cache_db}) : (),
	       $self->{special_protein_cache_dbhost} ? ("--cache-host", $self->{special_protein_cache_dbhost}) : (),
	       $self->{special_protein_cache_dbuser} ? ("--cache-user", $self->{special_protein_cache_dbuser}) : (),
	       $self->{special_protein_cache_dbpass} ? ("--cache-pass", $self->{special_protein_cache_dbpass}) : ());
    
    $ctx->stderr->log_cmd(@cmd);
    my $ok = run(\@cmd, '>', $out, $ctx->stderr->redirect);
    close($out);
    $results = [];
    if ($ok)
    {
	if (open(my $fh, "<", $out))
	{
	    my $l = <$fh>;
	    while (defined($l = <$fh>))
	    {
		chomp $l;
		push(@$results, [ split(/\t/, $l) ])
		}
	    close($fh);
	}
    }
    else
    {
	die "Error running rast_compute_specialty_genes ($?): @cmd\n" . $ctx->stderr->text_value;
    }
    
    #END compute_special_proteins
    my @_bad_returns;
    (ref($results) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"results\" (value was \"$results\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to compute_special_proteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($results);
}


=head2 annotate_special_proteins

  $genome_out = $obj->annotate_special_proteins($genome_in)

=over 4




=item Description


=back

=cut

sub annotate_special_proteins
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_special_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_special_proteins
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $dir = $self->{special_protein_dbdir};
    
    my $prots = File::Temp->new();
    close($prots);
    $genome_in->write_protein_translations_to_file($prots);
    
    my $out = File::Temp->new();
    
    my @opts;
    push(@opts, "--parallel", $self->{special_protein_threads});
    if ($self->{special_protein_cache_db})
    {
	push(@opts, "--cache-db", $self->{special_protein_cache_db});
	push(@opts,
	     $self->{special_protein_cache_dbhost} ? ("--cache-host", $self->{special_protein_cache_dbhost}) : (),
	     $self->{special_protein_cache_dbuser} ? ("--cache-user", $self->{special_protein_cache_dbuser}) : (),
	     $self->{special_protein_cache_dbpass} ? ("--cache-pass", $self->{special_protein_cache_dbpass}) : ());
    }
    else
    {
	push(@opts, "--no-cache");
    }
    
    my @cmd = ('rast_compute_specialty_genes',
	       "--in", $prots,
	       "--db-dir", $dir,
	       @opts);
    
    $ctx->stderr->log_cmd(@cmd);
    my $ok = run(\@cmd, '>', $out, $ctx->stderr->redirect);
    close($out);
    if ($ok)
    {
	if (open(my $fh, "<", $out))
	{
	    my $l = <$fh>;
	    while (defined($l = <$fh>))
	    {
		chomp $l;
		my($qid, $db, $sid, $qcov, $scov, $iden, $pvalue) = split(/\t/, $l);
		
		my $f = $genome_in->find_feature($qid);
		if (ref($f))
		{
		    push(@{$f->{similarity_associations}}, [$db, $sid, $qcov, $scov, $iden, $pvalue]);
		}
	    }
	    close($fh);
	}
    }
    else
    {
	die "Error running rast_compute_specialty_genes ($?): @cmd\n" . $ctx->stderr->text_value;
    }
    
    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    #END annotate_special_proteins
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_special_proteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 annotate_families_figfam_v1

  $genome_out = $obj->annotate_families_figfam_v1($genome_in)

=over 4




=item Description


=back

=cut

sub annotate_families_figfam_v1
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_families_figfam_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_families_figfam_v1
    
    #
    # Use the patric_call_proteins code from SEED to call figfamv1 families.
    # The path to this tool needs to be configured in the deployment configuration using
    # key patric_call_proteins_path.
    #
    
    my $remote = $self->{patric_call_proteins_remote_host};
    my $remote_user = $self->{patric_call_proteins_remote_user};
    my $exe = $self->{patric_call_proteins_path};
    
    if (!$remote && ($exe eq '' || ! -x $exe))
    {
	warn "PATRIC protein caller path not defined or invalid";
	$genome_out = $genome_in;
	goto do_return;
    }
    
    my $ff = $self->{patric_call_proteins_ff_path};
    my $md5_to_fam = $self->{patric_call_proteins_md5_to_fam_path};
    
    if (!$remote)
    {
	-d $ff or die "PATRIC protein caller figfam path not found";
	-f $md5_to_fam or die "PATRIC protein caller md5_to_fam path not found";
    }
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $prots = File::Temp->new();
    close($prots);
    $genome_in->write_protein_translations_to_file($prots);
    
    my $out = File::Temp->new();
    
    my(@cmd, $ok);
    if ($remote)
    {
	my $key = $self->{patric_call_proteins_remote_key};
	my @kopt = defined($key) ? ("-i", $key) : ();
	my @uopt = defined($remote_user) ? ("-l", $remote_user) : ();
	
	@cmd = ("ssh", @uopt, @kopt, $remote, "$exe -ff $ff --md5-to-fam $md5_to_fam -");
	print Dumper(\@cmd, $prots, $out);
	$ok = run(\@cmd, "<", "$prots", ">", $out, $ctx->stderr->redirect);
	close($out);
    }
    else
    {
	@cmd = ($exe, "-ff", $ff, "--md5-to-fam", $md5_to_fam, $prots);
	$ok = run(\@cmd, ">", $out, $ctx->stderr->redirect);
	close($out);
    }
    if (!$ok)
    {
	die "Error $? running @cmd\n" . $ctx->stderr->text_value;
    }
    
    open(my $fh, "<", $out) or die "Cannot open output temp $out: $!";
    
    while (<$fh>)
    {
	chomp;
	my($peg, $kfam, $kfun, @kscore) = split(/\t/);
	next unless $kfam;
	my $f = $genome_in->find_feature($peg);
	if (ref($f))
	{
	    push(@{$f->{family_assignments}}, [ 'FIGFAM', $kfam, $kfun ]);
	}
    }
    
    close($fh);
    
    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    do_return:
    #END annotate_families_figfam_v1
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_families_figfam_v1:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 annotate_families_patric

  $genome_out = $obj->annotate_families_patric($genome_in)

=over 4




=item Description


=back

=cut

sub annotate_families_patric
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_families_patric:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_families_patric
    
    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));
    
    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);
    
    my @cmd = ("place_proteins_into_pattyfams", "--output", $tmp_out,
	       $self->{patric_annotate_families_kmers},
	       $self->{patric_annotate_families_url},
	       $tmp_in);
    
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling patric families: $rc\non command @cmd";
    }
    
    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    
    #END annotate_families_patric
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_families_patric:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 annotate_families_patric_viral

  $genome_out = $obj->annotate_families_patric_viral($genome_in)

=over 4




=item Description


=back

=cut

sub annotate_families_patric_viral
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_families_patric_viral:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_families_patric_viral

    $genome_in = GenomeTypeObject->initialize($genome_in);

    #
    # Determine the genus for the given genome. We use that to find the appropriate PLF.
    #

    my $api = P3DataAPI->new();
    my $taxon_id = $genome_in->taxonomy_id();
    my @res = $api->query('taxonomy', ['eq', 'taxon_id', $taxon_id], ['select', 'lineage_ids', 'lineage_ranks']);

    my $genus;
    if (@res != 1)
    {
	warn "Could not find genus for $taxon_id\n";
    }
    else
    {
	my $ids = $res[0]->{lineage_ids};
	my $ranks = $res[0]->{lineage_ranks};
	while (@$ids)
	{
	    if ($ranks->[0] eq 'genus')
	    {
		$genus = $ids->[0];
		last;
	    }
	    shift @$ids;
	    shift @$ranks;
	}
    }
    
    #
    # Iterate over features and determine family IDs if present
    #
    my $plfam = $self->{viral_plfam}->{$genus};
    if ($genus && $plfam)
    {
	for my $f ($genome_in->features())
	{
	    if (my $ent = $plfam->{$f->{function}})
	    {
		my($type, $id, $product) = @$ent;
		push(@{$f->{family_assignments}}, [uc($type), $id, $product]);
	    }
	    elsif (my $ent = $self->{viral_pgfam}->{$f->{function}})
	    {
		my($type, $id, $product) = @$ent;
		push(@{$f->{family_assignments}}, [uc($type), $id, $product]);
	    }
	}
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
	    
    #END annotate_families_patric_viral
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_families_patric_viral:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 annotate_null_to_hypothetical

  $genome_out = $obj->annotate_null_to_hypothetical($genome_in)

=over 4




=item Description


=back

=cut

sub annotate_null_to_hypothetical
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_null_to_hypothetical:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_null_to_hypothetical
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $event = {
	tool_name => "annotate_null_to_hypothetical",
	execution_time => scalar gettimeofday,
	hostname => $self->{hostname},
    };
    my $event_id = $genome_in->add_analysis_event($event);
    
    for my $feature ($genome_in->features)
    {
	if ($feature->{function} eq '')
	{
	    $genome_in->update_function($ctx->user_id, $feature->{id}, "hypothetical protein", $event_id);
	}
    }
    
    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    #END annotate_null_to_hypothetical
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_null_to_hypothetical:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 remove_genbank_features

  $genome_out = $obj->remove_genbank_features($genome_in)

=over 4




=item Description


=back

=cut

sub remove_genbank_features
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to remove_genbank_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN remove_genbank_features

    #
    # Remove any features with a genbank_feature field. Used in pipelines
    # where we wish to import from genbank and propagate attributes, but not
    # use the genbank calls themselves.
    #

    $genome_in = GenomeTypeObject->initialize($genome_in);

    my $event = {
	tool_name => "remove_genbank_features",
	execution_time => scalar gettimeofday,
	hostname => $self->{hostname},
    };
    my $event_id = $genome_in->add_analysis_event($event);

    my $features = $genome_in->features;

    @$features = grep { ! $_->{genbank_feature} } @$features;
    
    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    #END remove_genbank_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to remove_genbank_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 annotate_strain_type_MLST

  $genome_out = $obj->annotate_strain_type_MLST($genome_in)

=over 4




=item Description


=back

=cut

sub annotate_strain_type_MLST
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_strain_type_MLST:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_strain_type_MLST
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    if (!defined($genome_in->{scientific_name}))
    {
	warn "No scientific name defined in annotate_strain_type_MLST";
    }
    else
    {
	my $contigs = $genome_in->extract_contig_sequences_to_temp_file();
	my $tmp_out = File::Temp->new();
	
	my $event = {
	    tool_name => "assign_st_to_genome",
	    execution_time => scalar gettimeofday,
	    hostname => $self->{hostname},
	};
	my $event_id = $genome_in->add_analysis_event($event);
	
	my @cmd = ("assign_st_to_genome",
		   "-s", 1, 
		   "-m", $self->{patric_mlst_dbdir},
		   "--org", $genome_in->{scientific_name});
	$ctx->stderr->log_cmd(@cmd);
	my $ok = run(\@cmd,
		     "<", $contigs,
		     ">", $tmp_out,
		     $ctx->stderr->redirect);
	
	close($tmp_out);
	unlink($contigs);
	
	if (open(my $fh, "<", $tmp_out))
	{
	    while (<$fh>)
	    {
		chomp;
		my($org, $db, $sids, $type, $coords, $loci, $profile, $tag) = split(/\t/);
		push(@{$genome_in->{typing}},
		 {
		     typing_method => 'MLST',
		     database => $db,
		     tag => $tag,
		     event_id => $event_id,
		 });
	    }
	}
    }	
    $genome_out = $genome_in->prepare_for_return();
    
    #END annotate_strain_type_MLST
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_strain_type_MLST:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 compute_cdd

  $return = $obj->compute_cdd($genome_in)

=over 4




=item Description


=back

=cut

sub compute_cdd
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to compute_cdd:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN compute_cdd
    #END compute_cdd
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to compute_cdd:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 annotate_proteins

  $return = $obj->annotate_proteins($genomeTO)

=over 4




=item Description


=back

=cut

sub annotate_proteins
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_proteins
    $return = $self->assign_functions_to_CDSs($genomeTO);
    #END annotate_proteins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 estimate_crude_phylogenetic_position_kmer

  $position_estimate = $obj->estimate_crude_phylogenetic_position_kmer($genomeTO)

=over 4




=item Description

Determine close genomes.
=back

=cut

sub estimate_crude_phylogenetic_position_kmer
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to estimate_crude_phylogenetic_position_kmer:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($position_estimate);
    #BEGIN estimate_crude_phylogenetic_position_kmer
    die "Not implemented";
    #END estimate_crude_phylogenetic_position_kmer
    my @_bad_returns;
    (!ref($position_estimate)) or push(@_bad_returns, "Invalid type for return variable \"position_estimate\" (value was \"$position_estimate\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to estimate_crude_phylogenetic_position_kmer:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($position_estimate);
}


=head2 find_close_neighbors

  $return = $obj->find_close_neighbors($genomeTO)

=over 4




=item Description


=back

=cut

sub find_close_neighbors
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to find_close_neighbors:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN find_close_neighbors

    my $kmer_service = Bio::KBase::KmerAnnotationByFigfam::Client->new($self->{kmer_service_url});

    if (! $kmer_service->can("estimate_closest_genomes"))
    {
	die "Kmer service does not support estimate_closest_genomes";
    }

    my @prots = map { [$_->{id}, $_->{function}, $_->{protein_translation}] }
    		grep { $_->{protein_translation} }
    		@{$genomeTO->{features}};
    
    my $dataset_name = $kmer_service->get_default_dataset_name();
    my $close = $kmer_service->estimate_closest_genomes(\@prots, undef);

    my $clist = $genomeTO->{close_genomes};
    if (!ref($clist))
    {
	$clist = [];
	$genomeTO->{close_genomes} = $clist;
    }

    push @$clist, { genome_id => $_->[0],
		    closeness_measure => $_->[1],
		    genome_name => $_->[2],
		    analysis_method => "KmerAnnotationByFigfam::estimate_closest_genomes:$dataset_name" } foreach @$close;

    $return = $genomeTO;

    #END find_close_neighbors
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to find_close_neighbors:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_strep_suis_repeat

  $return = $obj->call_features_strep_suis_repeat($genomeTO)

=over 4




=item Description

Interface to Strep repeats and "boxes" tools
=back

=cut

sub call_features_strep_suis_repeat
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_strep_suis_repeat:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_strep_suis_repeat

    $return = $self->_call_using_strep_repeats($genomeTO, "suis_repeat_annotation");

    #END call_features_strep_suis_repeat
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_strep_suis_repeat:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_strep_pneumo_repeat

  $return = $obj->call_features_strep_pneumo_repeat($genomeTO)

=over 4




=item Description


=back

=cut

sub call_features_strep_pneumo_repeat
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_strep_pneumo_repeat:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_strep_pneumo_repeat

    $return = $self->_call_using_strep_repeats($genomeTO, "pneumococcal_repeat_annotation");

    #END call_features_strep_pneumo_repeat
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_strep_pneumo_repeat:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 call_features_crispr

  $genome_out = $obj->call_features_crispr($genome_in)

=over 4




=item Description


=back

=cut

sub call_features_crispr
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_crispr

    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_crisprs", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id});
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling rRNAs: $rc\non command @cmd";
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));

    #END call_features_crispr
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 update_functions

  $genome_out = $obj->update_functions($genome_in, $functions, $event)

=over 4




=item Description


=back

=cut

sub update_functions
{
    my $self = shift;
    my($genome_in, $functions, $event) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($functions) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"functions\" (value was \"$functions\")");
    (ref($event) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"event\" (value was \"$event\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to update_functions:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN update_functions

    $genome_in = GenomeTypeObject->initialize($genome_in);

    #
    # If we don't have a valid event, create one here.
    #
    if (!$event->{tool_name})
    {
	$event = {
	    tool_name => "GenomeAnnotation::update_functions",
	    execution_time => scalar gettimeofday,
	    parameters => [],
	    hostname => $self->{hostname},
	};
    }
    my $event_id = $genome_in->add_analysis_event($event);

    for my $ent (@$functions)
    {
	my($fid, $function) = @$ent;
	$genome_in->update_function($ctx->user_id, $fid, $function, $event_id);
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();

    #END update_functions
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to update_functions:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 renumber_features

  $genome_out = $obj->renumber_features($genome_in)

=over 4




=item Description

Renumber features such that their identifiers are contiguous along contigs.

=back

=cut

sub renumber_features
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to renumber_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN renumber_features

    $genome_in = GenomeTypeObject->initialize($genome_in);

    my $event = {
	tool_name => "GenomeAnnotation::renumber_features",
	execution_time => scalar gettimeofday,
	parameters => [],
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);

    $genome_in->renumber_features($ctx->user_id, $event_id);

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();

    #END renumber_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to renumber_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 classify_amr

  $return = $obj->classify_amr($genomeTO)

=over 4




=item Description

Perform AMR classification.
=back

=cut

sub classify_amr
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to classify_amr:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN classify_amr

    my $genome_in = GenomeTypeObject->initialize($genomeTO);

    my $classifier = AdaboostClassify->new();

    my($genus, $species) = $genome_in->{scientific_name} =~ /^(\S+)\s+(\S+)/;
    my $name = "$genus $species";
    my $tmp = File::Temp->new();
    $tmp->close();
    $genome_in->write_contigs_to_file($tmp);
    my $res = $classifier->classify($name, "$tmp");

    my $idc = IDclient->new($genome_in);
    my $event = {
	tool_name => "AdaboostClassify",
	execution_time => scalar gettimeofday,
	parameters => [],
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);

    my $type = 'classifier_predicted_region';

    for my $classification (@$res)
    {
	my @flist;
	for my $f (@{$classification->{features}})
	{
	    my($contig, $start, $stop, $alpha, $round, $classifier, $function) = @$f;

	    my $len = $stop - $start + 1;
	    my $loc = [[$contig, $start, "+", $len]];
	    my $feat = $genome_in->add_feature({
		-id_client 	     => $idc,
		-id_prefix 	     => $genome_in->{id},
		-type 	     => $type,
		-location 	     => $loc,
		-function 	     => $function,
		-annotator => "$classifier classifier",
		-annotation      => "Classification by $classifier classifier with alpha=$alpha round=$round",
		-analysis_event_id 	     => $event_id,
	    });
	    push(@flist, [$feat->{id}, $alpha, $round, $function]);
	}
	my $cobj = {
	    name => $classification->{classifier},
	    comment => $classification->{comment},
	    antibiotics =>  [ split(/,\s+/, $classification->{antibiotic}) ],
	    accuracy => $classification->{accuracy},
	    area_under_roc_curve => $classification->{area_under_roc_curve},
	    f1_score => $classification->{f1_score},
	    cumulative_adaboost_score => $classification->{cumulative_adaboost_score},
	    sources => $classification->{sources},
	    sensitivity => $classification->{sensitivity},
	    event_id => $event_id,
	    features => [@flist],
	};
	push(@{$genome_in->{classifications}}, $cobj);
    }

    for my $f (<$tmp*>)
    {
	unlink($f);
    }

    $return = $genome_in->prepare_for_return();
    
    #END classify_amr
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to classify_amr:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 evaluate_genome

  $genome_out = $obj->evaluate_genome($genome_in, $params)

=over 4




=item Description

Perform genome evaluation.
=back

=cut

sub evaluate_genome
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to evaluate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN evaluate_genome

    #
    # Write to file before initializing.
    #

    my $tmpdir = File::Temp->newdir(undef, CLEANUP => 1);


    my $file = "$tmpdir/genome.gto";
    my $html = "GenomeReport.html";
    my $stdout;
    my $details = "genome_quality_details.txt";
    SeedUtils::write_encoded_object($genome_in, $file);

    my @ref;
    if (my $r = $params->{reference_genome_id})
    {
	@ref = ("--ref", $r);
    }
    #
    # Use our included path to BinningReports to find the webdetails file
    #
    my $rpt = $INC{'BinningReports.pm'};
    my $detail_template = dirname($rpt) . "/BinningReports/webdetails.tt";

    #
    # Check for new-style data
    #
    my @eval;
    my $eval_version;
    if (-d $self->{genome_evaluation_data})
    {
	@eval = ("--eval", $self->{genome_evaluation_data});
	if (open(my $fh, "<", "$self->{genome_evaluation_data}/VERSION"))
	{
	    $eval_version = <$fh>;
	    chomp $eval_version;
	    close($fh);
	    print "Evaluating genome with evaluation data version $eval_version\n";
	}
    }
    else
    {
	@eval = ("--predictors", $self->{genome_evaluation_predictors},
		 "--checkDir", $self->{genome_evaluation_checkg},
		 );
    }
    
    my @cmd = ("p3x-eval-gto",
	       @ref,
	       @eval,
	       "--deep",
	       "--template", $detail_template,
	       $file, $details, $html);

    print Dumper(\@cmd);
    
    my $ok = run(\@cmd, '>', \$stdout);

    if (!$ok)
    {
	die "Error $? running @cmd\n";
    }

    $genome_out = decode_json(scalar read_file("$file"));

    #END evaluate_genome
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to evaluate_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 export_genome

  $exported_data = $obj->export_genome($genome_in, $format, $feature_types)

=over 4




=item Description

Export genome typed object to one of the supported output formats:
genbank, embl, or gff.
If feature_types is a non-empty list, limit the output to the given
feature types.
=back

=cut

sub export_genome
{
    my $self = shift;
    my($genome_in, $format, $feature_types) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (!ref($format)) or push(@_bad_arguments, "Invalid type for argument \"format\" (value was \"$format\")");
    (ref($feature_types) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"feature_types\" (value was \"$feature_types\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to export_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($exported_data);
    #BEGIN export_genome

    my $coder = _get_coder();
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @type_flag = map { ("--feature-type", $_) } @$feature_types;

    if ($self->{special_protein_lookup_db})
    {
	push(@type_flag, "--specialty-gene-lookup-db", $self->{special_protein_lookup_db});
    }

    my @cmd = ("rast_export_genome", @type_flag, "--input", $tmp_in, "--output", $tmp_out, $format);

    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error exporting genome: $rc\non command @cmd";
    }

    $exported_data = scalar read_file("" . $tmp_out);

    #END export_genome
    my @_bad_returns;
    (!ref($exported_data)) or push(@_bad_returns, "Invalid type for return variable \"exported_data\" (value was \"$exported_data\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to export_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($exported_data);
}


=head2 enumerate_classifiers

  $return = $obj->enumerate_classifiers()

=over 4




=item Description

Enumerate the available classifiers. Returns the list of identifiers for
the classifiers.
=back

=cut

sub enumerate_classifiers
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN enumerate_classifiers

    $return = [];
    my $dir = $self->{kmer_classifier_data_directory};
    for my $ent (<$dir/*/groups>)
    {
	print STDERR "Try $ent\n";
	my($name) = $ent =~ m,([^/]+)/groups$,;
	push(@$return, $name);
    }
    #END enumerate_classifiers
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to enumerate_classifiers:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 query_classifier_groups

  $return = $obj->query_classifier_groups($classifier)

=over 4




=item Description

Query the groups included in the given classifier. This is a
mapping from the group name to the list of genome IDs included
in the group. Note that these are genome IDs native to the
system that created the classifier; currently these are
SEED genome IDs that may be translated using the
source IDs on the Genome entity.
=back

=cut

sub query_classifier_groups
{
    my $self = shift;
    my($classifier) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to query_classifier_groups:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN query_classifier_groups

    die "unimplemented";

    # disable # my $cobj = KmerClassifier->new("$self->{kmer_classifier_data_directory}/$classifier");
    # disable # $return = $cobj->group_membership_hash();

    #END query_classifier_groups
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to query_classifier_groups:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 query_classifier_taxonomies

  $return = $obj->query_classifier_taxonomies($classifier)

=over 4




=item Description

Query the taxonomy strings that this classifier maps.
=back

=cut

sub query_classifier_taxonomies
{
    my $self = shift;
    my($classifier) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to query_classifier_taxonomies:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN query_classifier_taxonomies
    #END query_classifier_taxonomies
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to query_classifier_taxonomies:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 classify_into_bins

  $return = $obj->classify_into_bins($classifier, $dna_input)

=over 4




=item Description

Classify a dataset, returning only the binned output.
=back

=cut

sub classify_into_bins
{
    my $self = shift;
    my($classifier, $dna_input) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    (ref($dna_input) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"dna_input\" (value was \"$dna_input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to classify_into_bins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN classify_into_bins

    return "unimplemented";
    
    # disable # my $cobj = KmerClassifier->new("$self->{kmer_classifier_data_directory}/$classifier");

    # disable # my $tmp = File::Temp->new;
    # disable # close($tmp);
    # disable # open(TMP, ">", $tmp) or die "cannot write tempfile $tmp: $!";
    # disable # for my $ent (@$dna_input)
# disable # {
    # disable # gjoseqlib::print_alignment_as_fasta(\*TMP, [$ent->[0], undef, $ent->[1]]);
# disable # }
    # disable # close(TMP);
    # disable # my($bins, $missed) = $cobj->classify("" . $tmp);
    # disable # 
    # disable # $return = $bins;
    
    #END classify_into_bins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to classify_into_bins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 classify_full

  $return_1, $raw_output, $unassigned = $obj->classify_full($classifier, $dna_input)

=over 4




=item Description

Classify a dataset, returning the binned output along with the raw assignments and the list of
sequences that were not assigned.
=back

=cut

sub classify_full
{
    my $self = shift;
    my($classifier, $dna_input) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    (ref($dna_input) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"dna_input\" (value was \"$dna_input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to classify_full:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return_1, $raw_output, $unassigned);
    #BEGIN classify_full


    die "unimplemented";
    # disable # my $cobj = KmerClassifier->new("$self->{kmer_classifier_data_directory}/$classifier");

    # disable # my $tmp = File::Temp->new;
    # disable # close($tmp);
    # disable # open(TMP, ">", $tmp) or die "cannot write tempfile $tmp: $!";
    
    # disable # my $raw_tmp = File::Temp->new();
    # disable # for my $ent (@$dna_input)
# disable # {
    # disable # gjoseqlib::print_alignment_as_fasta(\*TMP, [$ent->[0], undef, $ent->[1]]);
# disable # }
    # disable # close(TMP);
    # disable # my($bins, $missed) = $cobj->classify("" . $tmp, $raw_tmp);
    # disable # close($raw_tmp);
    
    # disable # $return_1 = $bins;
    # disable # $unassigned = $missed;
    # disable # $raw_output = read_file("" . $raw_tmp);
    
    #END classify_full
    my @_bad_returns;
    (ref($return_1) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return_1\" (value was \"$return_1\")");
    (!ref($raw_output)) or push(@_bad_returns, "Invalid type for return variable \"raw_output\" (value was \"$raw_output\")");
    (ref($unassigned) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"unassigned\" (value was \"$unassigned\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to classify_full:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return_1, $raw_output, $unassigned);
}


=head2 project_subsystems

  $genome_out = $obj->project_subsystems($genome_in)

=over 4




=item Description

Project subsystems.
=back

=cut

sub project_subsystems
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to project_subsystems:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN project_subsystems

    my $genome_in = GenomeTypeObject->initialize($genome_in);

    my $projector = $self->_subsystem_projector();

    my $parameters = [];
    my @keys = qw(subsystem-roles subsystem-variants subsystem-reference-data subsystem-variant-map);
    for my $k (@keys)
    {
	(my $sk = $k) =~ s/-/_/g;
	push(@$parameters, $sk => $self->{$sk});
    }

    my $event = {
	tool_name => "SubsystemProjector",
	execution_time => scalar gettimeofday,
	parameters => $parameters,
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);

    if ($projector)
    {
	my $ss = $projector->project_subsystems({genome_object => $genome_in, event_id => $event_id });
	$genome_in->{subsystems} = $ss;
    }
    else
    {
	warn "No subsystem projection parameters defined, skipping\n";
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();
    
    #END project_subsystems
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to project_subsystems:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 compute_genome_quality_control

  $genome_out = $obj->compute_genome_quality_control($genome_in)

=over 4




=item Description

Compute genome Quality Control scoring.
=back

=cut

sub compute_genome_quality_control
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to compute_genome_quality_control:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN compute_genome_quality_control

    my $coder = _get_coder();
    my $tmp_in = File::Temp->new(UNLINK => 1);
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new(UNLINK => 1);
    close($tmp_out);

    my @params = ("--genomeobj-file", "$tmp_in", "--new-genomeobj-file", "$tmp_out");
    my @cmd = ("rast2QC", @params);

    $ctx->stderr->log_cmd(@cmd);

    print STDERR "RUNNING @cmd\n";
    my $ok = run(\@cmd, $ctx->stderr->redirect);

    if (!$ok)
    {
	die "error running rast2QC: $?\non command @cmd\n" . $ctx->stderr->text_value;
    }

    $genome_out = GenomeTypeObject->new({file => "$tmp_out"});
    my $event = {
	tool_name => "rast2QC",
	execution_time => scalar gettimeofday,
	parameters => \@params,
	hostname => $self->{hostname},
    };

    my $event_id = $genome_out->add_analysis_event($event);

    $genome_out = $genome_out->prepare_for_return();

    #END compute_genome_quality_control
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to compute_genome_quality_control:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 default_workflow

  $return = $obj->default_workflow()

=over 4




=item Description


=back

=cut

sub default_workflow
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN default_workflow

    my @stages = (
	      { name => 'call_features_rRNA_SEED' },
	      { name => 'call_features_tRNA_trnascan' },
	      { name => 'call_features_repeat_region_SEED',
		    repeat_region_SEED_parameters => { } },
	      { name => 'call_selenoproteins', failure_is_not_fatal => 1 },
	      { name => 'call_pyrrolysoproteins', failure_is_not_fatal => 1 },
	      { name => 'call_features_strep_suis_repeat',
		condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
	      { name => 'call_features_strep_pneumo_repeat',
		condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
	      { name => 'call_features_crispr', failure_is_not_fatal => 1 },
	      { name => 'call_features_CDS_prodigal' },
	      { name => 'call_features_CDS_glimmer3', failure_is_not_fatal => 1, glimmer3_parameters => {} },
	      { name => 'prune_invalid_CDS_features', prune_invalid_CDS_features_parameters => {} },
	      { name => 'annotate_proteins_kmer_v2', kmer_v2_parameters => {} },
	      { name => 'annotate_proteins_kmer_v1', kmer_v1_parameters => { annotate_null_only => 1 }, failure_is_not_fatal => 1 },
	      { name => 'annotate_proteins_phage', phage_parameters => { annotate_null_only => 1 } },
	      { name => 'annotate_proteins_similarity', similarity_parameters => { annotate_null_only => 1 } },
	      { name => 'propagate_genbank_feature_metadata', propagate_genbank_feature_metadata_parameters => {} },
	      { name => 'resolve_overlapping_features', resolve_overlapping_features_parameters => {} },
	      { name => 'classify_amr', failure_is_not_fatal => 1 },
	      { name => 'renumber_features' },
              { name => 'annotate_special_proteins', failure_is_not_fatal => 1 },
	      { name => 'annotate_families_figfam_v1', failure_is_not_fatal => 1 },
	      { name => 'annotate_families_patric', failure_is_not_fatal => 1 },
	      { name => 'annotate_null_to_hypothetical' },
	      { name => 'find_close_neighbors', failure_is_not_fatal => 1 },
              { name => 'annotate_strain_type_MLST', failure_is_not_fatal => 1  },
	      # { name => 'call_features_prophage_phispy' },
		 );
    $return = { stages => \@stages };
    
    #END default_workflow
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to default_workflow:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 enumerate_recipes

  $recipes = $obj->enumerate_recipes()

=over 4




=item Description

Enumerate the loaded workflows. We always have a workflow named "default"; a
particular deployment of the genome annotation service may include additional workflows.
=back

=cut

sub enumerate_recipes
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($recipes);
    #BEGIN enumerate_recipes

    $recipes = [];

    my $def = $self->default_workflow();
    push(@$recipes, {
	id => 'default',
	name => 'Default annotation',
	description => "Standard PATRIC annotation",
	workflow => $def,
    });
    
    #
    # If we have a workflow directory configured, read from that.
    # Workflows are stored in workflow_dir/workflow-id/workflow.wf
    #
    my $dir = $self->{workflow_dir};
    if ($dir && opendir(my $dh, $dir))
    {
	my $coder = _get_coder();
	
	for my $id (sort { $a cmp $b } grep { -f "$dir/$_/workflow.wf" && $_ ne 'default' } readdir($dh))
	{
	    my $data;
	    eval {
		my $txt = read_file("$dir/$id/workflow.wf");
		$data = $coder->decode($txt);
	    };
	    if ($@)
	    {
		warn "error reading $dir/$id: $@";
	    }
	    else
	    {
		my $name = read_file("$dir/$id/name.txt");
		chomp $name;
		my $description = read_file("$dir/$id/description.txt");
		push(@$recipes, {
		    id => $id,
		    name => $name,
		    description => $description,
		    workflow => $data,
		});
	    }
	}
	closedir($dh);
    }
    
    #END enumerate_recipes
    my @_bad_returns;
    (ref($recipes) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"recipes\" (value was \"$recipes\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to enumerate_recipes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($recipes);
}


=head2 find_recipe

  $return = $obj->find_recipe($id)

=over 4




=item Description

Look up and return a particular named workflow.
=back

=cut

sub find_recipe
{
    my $self = shift;
    my($id) = @_;

    my @_bad_arguments;
    (!ref($id)) or push(@_bad_arguments, "Invalid type for argument \"id\" (value was \"$id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to find_recipe:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN find_recipe

    #
    # Support default workflow.
    #
    if ($id eq 'default')
    {
	my $def = $self->default_workflow();
	$return = {
	    id => 'default',
	    name => 'Default annotation',
	    description => "Standard PATRIC annotation",
	    workflow => $def,
	}
    }
    else
    {
	#
	# If we have a workflow directory configured, read from that.
	#
	my $dir = $self->{workflow_dir};
	#
	# Santize given workflow id to disallow loading from a path outside
	# the workflow directory
	#
	$id =~ s,/,,g;
	
	$return = {};
	
	if (open(my $fh, "<", "$dir/$id/workflow.wf"))
	{
	    my $coder = _get_coder();
	    
	    my $data;
	    eval {
		my $txt = read_file($fh, err_mode => 'quiet');
		$data = $coder->decode($txt);
	    };
	    if ($@)
	    {
		warn "error reading $dir/$id: $@";
	    }
	    else
	    {
		my $name = read_file("$dir/$id/name.txt", err_mode => 'quiet');
		chomp $name;
		my $description = read_file("$dir/$id/description.txt", err_mode => 'quiet');
		$return = {
		    id => $id,
		    name => $name,
		    description => $description,
		    workflow => $data,
		};
	    }
	    close($fh);
	}
    }

    #END find_recipe
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to find_recipe:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 complete_workflow_template

  $return = $obj->complete_workflow_template()

=over 4




=item Description

Return a workflow that includes all available stages. Not meant
(necessarily) for actual execution, but as a comprehensive list
of parts for users to use in assembling their own workflows.
=back

=cut

sub complete_workflow_template
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN complete_workflow_template

    my @stages = (
	      { name => 'call_features_rRNA_SEED' },
	      { name => 'call_features_tRNA_trnascan' },
	      { name => 'call_features_repeat_region_SEED',
		    repeat_region_SEED_parameters => { } },
	      { name => 'call_selenoproteins' },
	      { name => 'call_pyrrolysoproteins' },
	      { name => 'call_features_strep_suis_repeat',
		condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
	      { name => 'call_features_strep_pneumo_repeat',
		condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
	      { name => 'call_features_crispr', failure_is_not_fatal => 1 },
	      { name => 'call_features_CDS_glimmer3', glimmer3_parameters => { min_training_len => 2000 } },
	      { name => 'call_features_CDS_prodigal' },
	      { name => 'call_features_CDS_genemark' },
	      { name => 'annotate_proteins_kmer_v2', kmer_v2_parameters => { annotate_hypothetical_only => 0 } },
	      { name => 'annotate_proteins_kmer_v1', kmer_v1_parameters => { annotate_hypothetical_only => 1 } },
	      { name => 'annotate_proteins_similarity', similarity_parameters => { annotate_hypothetical_only => 1 } },
	      { name => 'annotate_null_to_hypothetical' },
	      { name => 'propagate_genbank_feature_metadata', propagate_genbank_feature_metadata_parameters => {} },
	      { name => 'resolve_overlapping_features', resolve_overlapping_features_parameters => {} },
	      { name => 'renumber_features' },
	      { name => 'find_close_neighbors', failure_is_not_fatal => 1 },
	      { name => 'call_features_prophage_phispy' },
		 );
    $return = { stages => \@stages };

    #END complete_workflow_template
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to complete_workflow_template:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 run_pipeline

  $genome_out = $obj->run_pipeline($genome_in, $workflow)

=over 4




=item Description


=back

=cut

sub run_pipeline
{
    my $self = shift;
    my($genome_in, $workflow) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($workflow) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"workflow\" (value was \"$workflow\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to run_pipeline:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN run_pipeline

    my %param_defs = (annotate_proteins_kmer_v1 => 'kmer_v1_parameters',
		      annotate_proteins_kmer_v2 => 'kmer_v2_parameters',
		      annotate_proteins_similarity => 'similarity_parameters',
		      annotate_proteins_phage => 'phage_parameters',
		      call_features_vigor4 => 'vigor4_parameters',
		      call_features_assembly_gap => 'assembly_gap_parameters',
		      call_features_repeat_region_SEED => 'repeat_regions_SEED_parameters',
		      call_features_CDS_glimmer3 => 'glimmer3_parameters',
		      call_features_repeat_region_SEED => 'repeat_region_SEED_parameters',
		      call_features_CDS_SEED_projection => 'SEED_projection_parameters',
		      call_features_ProtoCDS_kmer_v1 => 'kmer_v1_parameters',
		      prune_invalid_CDS_features => 'prune_invalid_CDS_features_parameters',
		      propagate_genbank_feature_metadata => 'propagate_genbank_feature_metadata_parameters',
		      split_gap_spanning_features => 'split_gap_spanning_features_params',
		      call_features_ProtoCDS_kmer_v2 => 'kmer_v2_parameters',
		      resolve_overlapping_features => 'resolve_overlapping_features_parameters',
		      propagate_genbank_feature_metadata => 'propagate_genbank_feature_metadata_parameters',
		      evaluate_genome => 'evaluate_genome_parameters',
		      );

    my $cur = $genome_in;
    for my $stage (@{$workflow->{stages}})
    {
	my $method = $stage->{name};
	my $condition = $stage->{condition};
	if ($condition)
	{
	    my $safe = Safe->new();
	    my $g = $safe->varglob('genome');
	    $$g = $cur;
	    my $ok = $safe->reval($condition);
	    print STDERR "Condition eval of '$condition' returns $ok\n";
	    if (!$ok)
	    {
		print STDERR "Skipping $method due to condition $condition\n";
		next;
	    }
	}
	my @params;
	if (my $param_def = $param_defs{$method})
	{
	    push(@params, $stage->{$param_def});
	}
	elsif ($method eq 'call_features_rRNA_SEED')
	{
	    # Special case.
	    push(@params, []);
	}
	if ($self->can($method))
	{
	    print STDERR "Call $method with " . Dumper(\@params);
	    print STDERR Dumper($stage);
	    my $out;
	    eval {
		$out = $self->$method($cur, @params);
	    };

	    if ($@)
	    {
		if ($stage->{failure_is_not_fatal})
		{
		    warn "Error invoking method $method: $@\nContinuing because failure_is_not_fatal flag is set";
		    if (ref($cur) && ref($cur) ne 'HASH')
		    {
			$cur = $cur->prepare_for_return;
		    }
		}
		else
		{
		    die "Error invoking method $method: $@";
		}
	    }
	    else
	    {
		print STDERR "Finished\n";
		# print Dumper($out);
		$cur = $out;
	    }
	}
	else
	{
	    die "Trying to call invalid method $method";
	}
    }

    $genome_out = $cur;

    #END run_pipeline
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to run_pipeline:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($genome_out);
}


=head2 pipeline_batch_start

  $batch_id = $obj->pipeline_batch_start($genomes, $workflow)

=over 4




=item Description


=back

=cut

sub pipeline_batch_start
{
    my $self = shift;
    my($genomes, $workflow) = @_;

    my @_bad_arguments;
    (ref($genomes) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"genomes\" (value was \"$genomes\")");
    (ref($workflow) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"workflow\" (value was \"$workflow\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to pipeline_batch_start:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($batch_id);
    #BEGIN pipeline_batch_start

    #
    # Construct an AWE workflow and submit it to our AWE server.
    #
    # The genomes have already been uploaded - the genomes list here is a
    # list of handle objects.
    #
    # The workflow we create will have one task per genome submitted. The data
    # passed to each task are the serialized form of the handle and the
    # workflow description.
    #

    my $json = JSON::XS->new->pretty(1);
    my $shock = Bio::KBase::GenomeAnnotation::Shock->new($self->{shock_server}, $ctx->token);
    my $awe = Bio::KBase::GenomeAnnotation::Awe->new($self->{awe_server}, $ctx->token);
    
    my $job = Bio::KBase::GenomeAnnotation::Awe::JobDescription->new(pipeline => 'rasttk',
								     name => 'rasttk',
								     project => 'rasttk',
								     user => $ctx->user_id,
								     clientgroups => '');
    
    my $i = 0;
    for my $gspec (@$genomes)
    {
	my($genome_id, $handle, $filename) = @$gspec{'genome_id', 'data', 'filename'};
	
	my $txt = $json->encode([$handle, $workflow]);
	my $node = $shock->put_file_data($txt);
	
	my $in_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("pipeinput_$i.json", $shock->server, $node);
	my $out_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("pipeoutput_$i.json", $shock->server);
	my $stdout_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("stdout_$i.txt", $shock->server);
	my $stderr_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("stderr_$i.txt", $shock->server);
	$i++;
	
	my $id = $job->add_task("rast_run_pipeline_local",
				"rast_run_pipeline_local",
				join(" ",
				     $in_file->in_name,
				     $out_file->name,
				     $stdout_file->name,
				     $stderr_file->name),
				[],
				[$in_file],
				[$out_file, $stdout_file, $stderr_file],
				undef, undef, $awe,
			        { genome_id => $genome_id, filename => $filename });
    }
    
    $batch_id = $awe->submit($job);
    print STDERR "submitted $batch_id\n";
    
    #END pipeline_batch_start
    my @_bad_returns;
    (!ref($batch_id)) or push(@_bad_returns, "Invalid type for return variable \"batch_id\" (value was \"$batch_id\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to pipeline_batch_start:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($batch_id);
}


=head2 pipeline_batch_status

  $status = $obj->pipeline_batch_status($batch_id)

=over 4




=item Description


=back

=cut

sub pipeline_batch_status
{
    my $self = shift;
    my($batch_id) = @_;

    my @_bad_arguments;
    (!ref($batch_id)) or push(@_bad_arguments, "Invalid type for argument \"batch_id\" (value was \"$batch_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to pipeline_batch_status:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($status);
    #BEGIN pipeline_batch_status

    my $awe = Bio::KBase::GenomeAnnotation::Awe->new($self->{awe_server}, $ctx->token);
    
    my $job = $awe->job($batch_id);

    #
    # We could check to see if the job's user matches the authenticated user here. I
    # don't think it's necessary at this point.
    #

    #
    # Loop over the tasks collecting the status / file information to be formatted and returned.
    #

    my $tasks = $job->{tasks};

    my $details = [];

    for my $task (@$tasks)
    {
	my $tinfo = {
	    genome_id => $task->{userattr}->{genome_id},
	    filename => $task->{userattr}->{filename},
	    status => $task->{state},
	    creation_date => $task->{createddate},
	    start_date => $task->{starteddate},
	    completion_date => $task->{completeddate},
	};

	my $outfiles = $task->{outputs};
	#
	# Support both old and new format for file list by constructing
	# a new format list from the old.
	#
	if (ref($outfiles) eq 'HASH')
	{
	    my $new = [];
	    while (my($file, $hash) = each %$outfiles)
	    {
		my $nhash = { filename => $file, %$hash };
		push(@$new, $nhash);
	    }
	    $outfiles = $new;
	}
	    
	for my $f (@$outfiles)
	{
	    my $key;
	    my $filename = $f->{filename};
	    if ($filename =~ /^stderr/) { $key = 'stderr'; }
	    elsif ($filename =~ /^stdout/ ) { $key = 'stdout'; }
	    elsif ($filename =~ /^pipeoutput/) { $key = 'output'; }
	    else { next; }
	
	    if ($f->{url} ne '')
	    {
		$tinfo->{$key} = {
		    file_name => $filename,
		    id => $f->{node},
		    type => 'shock',
		    url => $f->{host},
		    remote_md5 => '',
		    remote_sha1 => '',
		};
	    }
	}
	push (@$details, $tinfo);
    }

    $status = {
	status => $job->{info}->{state},
	submit_date => $job->{info}->{submittime},
	start_date => $job->{info}->{startedtime},
	completion_date => $job->{info}->{completedtime},
	details => $details,
    };
    
    #END pipeline_batch_status
    my @_bad_returns;
    (ref($status) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"status\" (value was \"$status\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to pipeline_batch_status:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($status);
}


=head2 pipeline_batch_enumerate_batches

  $batches = $obj->pipeline_batch_enumerate_batches()

=over 4




=item Description


=back

=cut

sub pipeline_batch_enumerate_batches
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($batches);
    #BEGIN pipeline_batch_enumerate_batches

    #
    # Call AWE to find all tasks with project=rasttk and user=logged-in user
    #
    
    my $awe = Bio::KBase::GenomeAnnotation::Awe->new($self->{awe_server}, $ctx->token);
    my($list, $error) = $awe->enumerate_jobs({ "info.project" => 'rasttk', "info.user" => $ctx->user_id });

    $batches = [];
    if (ref($list))
    {
	$batches = [ map { [ $_->{id}, $_->{info}->{submittime} ] } @$list ];
    }
    else
    {
	die "pipeline_batch_enumerate_batches task query error: $error";
    }

    #END pipeline_batch_enumerate_batches
    my @_bad_returns;
    (ref($batches) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"batches\" (value was \"$batches\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to pipeline_batch_enumerate_batches:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($batches);
}





=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}


1;
