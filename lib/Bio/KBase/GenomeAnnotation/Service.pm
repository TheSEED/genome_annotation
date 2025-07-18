package Bio::KBase::GenomeAnnotation::Service;


use strict;
use Data::Dumper;
use Moose;
use POSIX;
use JSON;
use File::Temp;
use File::Slurp;
use Class::Load qw();
use Config::Simple;
use Time::HiRes qw(gettimeofday);


use P3AuthToken;
use P3TokenValidator;

my $g_hostname = `hostname`;
chomp $g_hostname;
$g_hostname ||= 'unknown-host';

extends 'RPC::Any::Server::JSONRPC::PSGI';

has 'instance_dispatch' => (is => 'ro', isa => 'HashRef');
has 'user_auth' => (is => 'ro', isa => 'UserAuth');
has 'valid_methods' => (is => 'ro', isa => 'HashRef', lazy => 1,
			builder => '_build_valid_methods');
has 'validator' => (is => 'ro', isa => 'P3TokenValidator', lazy => 1, builder => '_build_validator');
our $CallContext;

our %return_counts = (
        'genome_ids_to_genomes' => 1,
        'create_genome' => 1,
        'create_genome_from_genbank' => 1,
        'create_genome_from_SEED' => 1,
        'create_genome_from_RAST' => 1,
        'set_metadata' => 1,
        'add_contigs' => 1,
        'add_contigs_from_handle' => 1,
        'import_sra_metadata' => 1,
        'compute_sars2_variation' => 1,
        'add_features' => 1,
        'genomeTO_to_reconstructionTO' => 1,
        'genomeTO_to_feature_data' => 1,
        'reconstructionTO_to_roles' => 1,
        'reconstructionTO_to_subsystems' => 1,
        'assign_functions_to_CDSs' => 1,
        'annotate_genome' => 1,
        'call_selenoproteins' => 1,
        'call_pyrrolysoproteins' => 1,
        'call_features_selenoprotein' => 1,
        'call_features_pyrrolysoprotein' => 1,
        'call_features_insertion_sequences' => 1,
        'call_features_rRNA_SEED' => 1,
        'call_features_tRNA_trnascan' => 1,
        'call_RNAs' => 1,
        'call_features_lowvan' => 1,
        'call_features_vigor4' => 1,
        'call_features_vipr_mat_peptide' => 1,
        'call_features_CDS_glimmer3' => 1,
        'call_features_CDS_prodigal' => 1,
        'call_features_CDS_genemark' => 1,
        'call_features_CDS_phanotate' => 1,
        'prune_invalid_CDS_features' => 1,
        'call_features_CDS_SEED_projection' => 1,
        'call_features_CDS_FragGeneScan' => 1,
        'call_features_repeat_region_SEED' => 1,
        'call_features_prophage_phispy' => 1,
        'call_features_scan_for_matches' => 1,
        'call_features_assembly_gap' => 1,
        'split_gap_spanning_features' => 1,
        'translate_untranslated_proteins' => 1,
        'annotate_proteins_similarity' => 1,
        'annotate_proteins_phage' => 1,
        'annotate_proteins_kmer_v1' => 1,
        'annotate_proteins_kmer_v2' => 1,
        'resolve_overlapping_features' => 1,
        'propagate_genbank_feature_metadata' => 1,
        'call_features_ProtoCDS_kmer_v1' => 1,
        'call_features_ProtoCDS_kmer_v2' => 1,
        'enumerate_special_protein_databases' => 1,
        'compute_special_proteins' => 1,
        'annotate_special_proteins' => 1,
        'annotate_families_figfam_v1' => 1,
        'annotate_families_patric' => 1,
        'annotate_families_patric_viral' => 1,
        'annotate_null_to_hypothetical' => 1,
        'remove_genbank_features' => 1,
        'annotate_strain_type_MLST' => 1,
        'compute_cdd' => 1,
        'annotate_proteins' => 1,
        'estimate_crude_phylogenetic_position_kmer' => 1,
        'find_close_neighbors' => 1,
        'call_features_strep_suis_repeat' => 1,
        'call_features_strep_pneumo_repeat' => 1,
        'call_features_crispr' => 1,
        'update_functions' => 1,
        'renumber_features' => 1,
        'classify_amr' => 1,
        'evaluate_genome' => 1,
        'export_genome' => 1,
        'enumerate_classifiers' => 1,
        'query_classifier_groups' => 1,
        'query_classifier_taxonomies' => 1,
        'classify_into_bins' => 1,
        'classify_full' => 3,
        'project_subsystems' => 1,
        'compute_genome_quality_control' => 1,
        'default_workflow' => 1,
        'enumerate_recipes' => 1,
        'find_recipe' => 1,
        'complete_workflow_template' => 1,
        'run_pipeline' => 1,
        'pipeline_batch_start' => 1,
        'pipeline_batch_status' => 1,
        'pipeline_batch_enumerate_batches' => 1,
        'version' => 1,
);

our %method_authentication = (
        'genome_ids_to_genomes' => 'none',
        'create_genome' => 'none',
        'create_genome_from_genbank' => 'none',
        'create_genome_from_SEED' => 'none',
        'create_genome_from_RAST' => 'none',
        'set_metadata' => 'none',
        'add_contigs' => 'none',
        'add_contigs_from_handle' => 'none',
        'import_sra_metadata' => 'none',
        'compute_sars2_variation' => 'none',
        'add_features' => 'none',
        'genomeTO_to_reconstructionTO' => 'none',
        'genomeTO_to_feature_data' => 'none',
        'reconstructionTO_to_roles' => 'none',
        'reconstructionTO_to_subsystems' => 'none',
        'assign_functions_to_CDSs' => 'none',
        'annotate_genome' => 'none',
        'call_selenoproteins' => 'none',
        'call_pyrrolysoproteins' => 'none',
        'call_features_selenoprotein' => 'none',
        'call_features_pyrrolysoprotein' => 'none',
        'call_features_insertion_sequences' => 'none',
        'call_features_rRNA_SEED' => 'none',
        'call_features_tRNA_trnascan' => 'none',
        'call_RNAs' => 'none',
        'call_features_lowvan' => 'none',
        'call_features_vigor4' => 'none',
        'call_features_vipr_mat_peptide' => 'none',
        'call_features_CDS_glimmer3' => 'none',
        'call_features_CDS_prodigal' => 'none',
        'call_features_CDS_genemark' => 'none',
        'call_features_CDS_phanotate' => 'none',
        'prune_invalid_CDS_features' => 'none',
        'call_features_CDS_SEED_projection' => 'none',
        'call_features_CDS_FragGeneScan' => 'none',
        'call_features_repeat_region_SEED' => 'none',
        'call_features_prophage_phispy' => 'none',
        'call_features_scan_for_matches' => 'none',
        'call_features_assembly_gap' => 'none',
        'split_gap_spanning_features' => 'none',
        'translate_untranslated_proteins' => 'none',
        'annotate_proteins_similarity' => 'none',
        'annotate_proteins_phage' => 'none',
        'annotate_proteins_kmer_v1' => 'none',
        'annotate_proteins_kmer_v2' => 'none',
        'resolve_overlapping_features' => 'none',
        'propagate_genbank_feature_metadata' => 'none',
        'call_features_ProtoCDS_kmer_v1' => 'none',
        'call_features_ProtoCDS_kmer_v2' => 'none',
        'enumerate_special_protein_databases' => 'none',
        'compute_special_proteins' => 'none',
        'annotate_special_proteins' => 'none',
        'annotate_families_figfam_v1' => 'none',
        'annotate_families_patric' => 'none',
        'annotate_families_patric_viral' => 'none',
        'annotate_null_to_hypothetical' => 'none',
        'remove_genbank_features' => 'none',
        'annotate_strain_type_MLST' => 'none',
        'compute_cdd' => 'none',
        'annotate_proteins' => 'none',
        'estimate_crude_phylogenetic_position_kmer' => 'none',
        'find_close_neighbors' => 'none',
        'call_features_strep_suis_repeat' => 'none',
        'call_features_strep_pneumo_repeat' => 'none',
        'call_features_crispr' => 'none',
        'update_functions' => 'none',
        'renumber_features' => 'none',
        'classify_amr' => 'none',
        'evaluate_genome' => 'none',
        'export_genome' => 'none',
        'enumerate_classifiers' => 'none',
        'query_classifier_groups' => 'none',
        'query_classifier_taxonomies' => 'none',
        'classify_into_bins' => 'none',
        'classify_full' => 'none',
        'project_subsystems' => 'none',
        'compute_genome_quality_control' => 'none',
        'default_workflow' => 'none',
        'enumerate_recipes' => 'none',
        'find_recipe' => 'none',
        'complete_workflow_template' => 'none',
        'run_pipeline' => 'none',
        'pipeline_batch_start' => 'required',
        'pipeline_batch_status' => 'required',
        'pipeline_batch_enumerate_batches' => 'required',
);

sub _build_validator
{
    my($self) = @_;
    return P3TokenValidator->new();

}


sub _build_valid_methods
{
    my($self) = @_;
    my $methods = {
        'genome_ids_to_genomes' => 1,
        'create_genome' => 1,
        'create_genome_from_genbank' => 1,
        'create_genome_from_SEED' => 1,
        'create_genome_from_RAST' => 1,
        'set_metadata' => 1,
        'add_contigs' => 1,
        'add_contigs_from_handle' => 1,
        'import_sra_metadata' => 1,
        'compute_sars2_variation' => 1,
        'add_features' => 1,
        'genomeTO_to_reconstructionTO' => 1,
        'genomeTO_to_feature_data' => 1,
        'reconstructionTO_to_roles' => 1,
        'reconstructionTO_to_subsystems' => 1,
        'assign_functions_to_CDSs' => 1,
        'annotate_genome' => 1,
        'call_selenoproteins' => 1,
        'call_pyrrolysoproteins' => 1,
        'call_features_selenoprotein' => 1,
        'call_features_pyrrolysoprotein' => 1,
        'call_features_insertion_sequences' => 1,
        'call_features_rRNA_SEED' => 1,
        'call_features_tRNA_trnascan' => 1,
        'call_RNAs' => 1,
        'call_features_lowvan' => 1,
        'call_features_vigor4' => 1,
        'call_features_vipr_mat_peptide' => 1,
        'call_features_CDS_glimmer3' => 1,
        'call_features_CDS_prodigal' => 1,
        'call_features_CDS_genemark' => 1,
        'call_features_CDS_phanotate' => 1,
        'prune_invalid_CDS_features' => 1,
        'call_features_CDS_SEED_projection' => 1,
        'call_features_CDS_FragGeneScan' => 1,
        'call_features_repeat_region_SEED' => 1,
        'call_features_prophage_phispy' => 1,
        'call_features_scan_for_matches' => 1,
        'call_features_assembly_gap' => 1,
        'split_gap_spanning_features' => 1,
        'translate_untranslated_proteins' => 1,
        'annotate_proteins_similarity' => 1,
        'annotate_proteins_phage' => 1,
        'annotate_proteins_kmer_v1' => 1,
        'annotate_proteins_kmer_v2' => 1,
        'resolve_overlapping_features' => 1,
        'propagate_genbank_feature_metadata' => 1,
        'call_features_ProtoCDS_kmer_v1' => 1,
        'call_features_ProtoCDS_kmer_v2' => 1,
        'enumerate_special_protein_databases' => 1,
        'compute_special_proteins' => 1,
        'annotate_special_proteins' => 1,
        'annotate_families_figfam_v1' => 1,
        'annotate_families_patric' => 1,
        'annotate_families_patric_viral' => 1,
        'annotate_null_to_hypothetical' => 1,
        'remove_genbank_features' => 1,
        'annotate_strain_type_MLST' => 1,
        'compute_cdd' => 1,
        'annotate_proteins' => 1,
        'estimate_crude_phylogenetic_position_kmer' => 1,
        'find_close_neighbors' => 1,
        'call_features_strep_suis_repeat' => 1,
        'call_features_strep_pneumo_repeat' => 1,
        'call_features_crispr' => 1,
        'update_functions' => 1,
        'renumber_features' => 1,
        'classify_amr' => 1,
        'evaluate_genome' => 1,
        'export_genome' => 1,
        'enumerate_classifiers' => 1,
        'query_classifier_groups' => 1,
        'query_classifier_taxonomies' => 1,
        'classify_into_bins' => 1,
        'classify_full' => 1,
        'project_subsystems' => 1,
        'compute_genome_quality_control' => 1,
        'default_workflow' => 1,
        'enumerate_recipes' => 1,
        'find_recipe' => 1,
        'complete_workflow_template' => 1,
        'run_pipeline' => 1,
        'pipeline_batch_start' => 1,
        'pipeline_batch_status' => 1,
        'pipeline_batch_enumerate_batches' => 1,
        'version' => 1,
    };
    return $methods;
}

#
# Override method from RPC::Any::Server::JSONRPC 
# to eliminate the deprecation warning for Class::MOP::load_class.
#
sub _default_error {
    my ($self, %params) = @_;
    my $version = $self->default_version;
    $version =~ s/\./_/g;
    my $error_class = "JSON::RPC::Common::Procedure::Return::Version_${version}::Error";
    Class::Load::load_class($error_class);
    my $error = $error_class->new(%params);
    my $return_class = "JSON::RPC::Common::Procedure::Return::Version_$version";
    Class::Load::load_class($return_class);
    return $return_class->new(error => $error);
}


#override of RPC::Any::Server
sub handle_error {
    my ($self, $error) = @_;
    
    unless (ref($error) eq 'HASH' ||
           (blessed $error and $error->isa('RPC::Any::Exception'))) {
        $error = RPC::Any::Exception::PerlError->new(message => $error);
    }
    my $output;
    eval {
        my $encoded_error = $self->encode_output_from_exception($error);
        $output = $self->produce_output($encoded_error);
    };
    
    return $output if $output;
    
    die "$error\n\nAlso, an error was encountered while trying to send"
        . " this error: $@\n";
}

#override of RPC::Any::JSONRPC
sub encode_output_from_exception {
    my ($self, $exception) = @_;
    my %error_params;
    if (ref($exception) eq 'HASH') {
        %error_params = %{$exception};
        if(defined($error_params{context})) {
            my @errlines;
            $errlines[0] = $error_params{message};
            push @errlines, split("\n", $error_params{data});
            delete $error_params{context};
        }
    } else {
        %error_params = (
            message => $exception->message,
            code    => $exception->code,
        );
    }
    my $json_error;
    if ($self->_last_call) {
        $json_error = $self->_last_call->return_error(%error_params);
    }
    # Default to default_version. This happens when we throw an exception
    # before inbound parsing is complete.
    else {
        $json_error = $self->_default_error(%error_params);
    }
    return $self->encode_output_from_object($json_error);
}

#
# another override.
#
sub get_package_isa {
    my ($self, $module) = @_;
    my $original_isa;
    { no strict 'refs'; $original_isa = \@{"${module}::ISA"}; }
    my @new_isa = @$original_isa;

    my $base = $self->package_base;
    if (not $module->isa($base)) {
        Class::Load::load_class($base);
        push(@new_isa, $base);
    }
    return \@new_isa;
}
sub trim {
    my ($str) = @_;
    if (!(defined $str)) {
        return $str;
    }
    $str =~ s/^\s+|\s+$//g;
    return $str;
}

sub getIPAddress {
    my ($self) = @_;
    my $xFF = trim($self->_plack_req->header("X-Forwarded-For"));
    my $realIP = trim($self->_plack_req->header("X-Real-IP"));
    # my $nh = $self->config->{"dont_trust_x_ip_headers"};
    my $nh;
    my $trustXHeaders = !(defined $nh) || $nh ne "true";

    if ($trustXHeaders) {
        if ($xFF) {
            my @tmp = split(",", $xFF);
            return trim($tmp[0]);
        }
        if ($realIP) {
            return $realIP;
        }
    }
    return $self->_plack_req->address;
}

#
# Ping method reflected from /ping on the service.
#
sub ping
{
    my($self, $env) = @_;
    return [ 200, ["Content-type" => "text/plain"], [ "OK\n" ] ];
}


#
# Authenticated ping method reflected from /auth_ping on the service.
#
sub auth_ping
{
    my($self, $env) = @_;

    my $req = Plack::Request->new($env);
    my $token = $req->header("Authorization");

    if (!$token)
    {
	return [401, [], ["Authentication required"]];
    }

    my $auth_token = P3AuthToken->new(token => $token, ignore_authrc => 1);
    my($valid, $validate_err) = $self->validator->validate($auth_token);

    if ($valid)
    {
	return [200, ["Content-type" => "text/plain"], ["OK " . $auth_token->user_id . "\n"]];
    }
    else
    {
        warn "Token validation error $validate_err\n";
	return [403, [], "Authentication failed"];
    }
}

sub call_method {
    my ($self, $data, $method_info) = @_;

    my ($module, $method, $modname) = @$method_info{qw(module method modname)};
    
    my $ctx = Bio::KBase::GenomeAnnotation::ServiceContext->new(client_ip => $self->getIPAddress());
    $ctx->module($modname);
    $ctx->method($method);
    $ctx->call_id($self->{_last_call}->{id});
    
    my $args = $data->{arguments};

{
    # Service GenomeAnnotation requires authentication.

    my $method_auth = $method_authentication{$method};
    $ctx->authenticated(0);
    if ($method_auth eq 'none')
    {
	# No authentication required here. Move along.
    }
    else
    {
	my $token = $self->_plack_req->header("Authorization");

	if (!$token && $method_auth eq 'required')
	{
	    $self->exception('PerlError', "Authentication required for GenomeAnnotation but no authentication header was passed");
	}

	my $auth_token = P3AuthToken->new(token => $token, ignore_authrc => 1);
	my($valid, $validate_err) = $self->validator->validate($auth_token);
	# Only throw an exception if authentication was required and it fails
	if ($method_auth eq 'required' && !$valid)
	{
	    $self->exception('PerlError', "Token validation failed: $validate_err");
	} elsif ($valid) {
	    $ctx->authenticated(1);
	    $ctx->user_id($auth_token->user_id);
	    $ctx->token( $token);
	}
    }
}
    my $new_isa = $self->get_package_isa($module);
    no strict 'refs';
    local @{"${module}::ISA"} = @$new_isa;
    local $CallContext = $ctx;
    my @result;
    {
	# 
	# Process tag and metadata information if present.
	#
	my $tag = $self->_plack_req->header("Kbrpc-Tag");
	if (!$tag)
	{
	    $self->{hostname} ||= $g_hostname;

	    my ($t, $us) = gettimeofday;
	    $us = sprintf("%06d", $us);
	    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	    $tag = "S:$self->{hostname}:$$:$ts";
	}
	local $ENV{KBRPC_TAG} = $tag;
	my $kb_metadata = $self->_plack_req->header("Kbrpc-Metadata");
	my $kb_errordest = $self->_plack_req->header("Kbrpc-Errordest");
	local $ENV{KBRPC_METADATA} = $kb_metadata if $kb_metadata;
	local $ENV{KBRPC_ERROR_DEST} = $kb_errordest if $kb_errordest;

	my $stderr = Bio::KBase::GenomeAnnotation::ServiceStderrWrapper->new($ctx);
	$ctx->stderr($stderr);

	#
	# Set up environment for user-level error reporting.
	#
	my $user_error = File::Temp->new(UNLINK => 1);
	close($user_error);
	$ENV{P3_USER_ERROR_DESTINATION} = "$user_error";

        my $xFF = $self->_plack_req->header("X-Forwarded-For");
	
        my $err;
        eval {
	    local $SIG{__WARN__} = sub {
		my($msg) = @_;
		print STDERR $msg;
	    };

            @result = $module->$method(@{ $data->{arguments} });
        };
	
        if ($@)
        {
            my $err = $@;
	    $stderr->log($err);
	    $ctx->stderr(undef);
	    undef $stderr;
            my $nicerr;
	    my $str = "$err";
	    my $msg = $str;
	    $msg =~ s/ at [^\s]+.pm line \d+.\n$//;
	    
	    # If user-level error present, replace message with that
	    if (-s "$user_error")
	    {
	        $msg = read_file("$user_error");
		$str = $msg;
	    }
	    $nicerr =  {code => -32603, # perl error from RPC::Any::Exception
                            message => $msg,
                            data => $str,
                            context => $ctx
                            };
            die $nicerr;
        }
	$ctx->stderr(undef);
	undef $stderr;
    }
    my $result;
    if ($return_counts{$method} == 1)
    {
        $result = [[$result[0]]];
    }
    else
    {
        $result = \@result;
    }
    return $result;
}


sub get_method
{
    my ($self, $data) = @_;
    
    my $full_name = $data->{method};
    
    $full_name =~ /^(\S+)\.([^\.]+)$/;
    my ($package, $method) = ($1, $2);
    
    if (!$package || !$method) {
	$self->exception('NoSuchMethod',
			 "'$full_name' is not a valid method. It must"
			 . " contain a package name, followed by a period,"
			 . " followed by a method name.");
    }

    if (!$self->valid_methods->{$method})
    {
	$self->exception('NoSuchMethod',
			 "'$method' is not a valid method in service GenomeAnnotation.");
    }
	
    my $inst = $self->instance_dispatch->{$package};
    my $module;
    if ($inst)
    {
	$module = $inst;
    }
    else
    {
	$module = $self->get_module($package);
	if (!$module) {
	    $self->exception('NoSuchMethod',
			     "There is no method package named '$package'.");
	}
	
	Class::Load::load_class($module);
    }
    
    if (!$module->can($method)) {
	$self->exception('NoSuchMethod',
			 "There is no method named '$method' in the"
			 . " '$package' package.");
    }
    
    return { module => $module, method => $method, modname => $package };
}

package Bio::KBase::GenomeAnnotation::ServiceContext;

use strict;

=head1 NAME

Bio::KBase::GenomeAnnotation::ServiceContext

head1 DESCRIPTION

A KB RPC context contains information about the invoker of this
service. If it is an authenticated service the authenticated user
record is available via $context->user. The client IP address
is available via $context->client_ip.

=cut

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(user_id client_ip authenticated token
                             module method call_id hostname stderr));

sub new
{
    my($class, @opts) = @_;

    if (!defined($opts[0]) || ref($opts[0]))
    {
        # We were invoked by old code that stuffed a logger in here.
	# Strip that option.
	shift @opts;
    }
    
    my $self = {
        hostname => $g_hostname,
        @opts,
    };

    return bless $self, $class;
}

package Bio::KBase::GenomeAnnotation::ServiceStderrWrapper;

use strict;
use POSIX;
use Time::HiRes 'gettimeofday';

sub new
{
    my($class, $ctx) = @_;
    my $self = {
    };
    my $dest = $ENV{KBRPC_ERROR_DEST} if exists $ENV{KBRPC_ERROR_DEST};
    my $tag = $ENV{KBRPC_TAG} if exists $ENV{KBRPC_TAG};
    my ($t, $us) = gettimeofday();
    $us = sprintf("%06d", $us);
    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);

    my $name = join(".", $ctx->module, $ctx->method, $ctx->hostname, $ts);

    if ($dest && $dest =~ m,^/,)
    {
	#
	# File destination
	#
	my $fh;

	if ($tag)
	{
	    $tag =~ s,/,_,g;
	    $dest = "$dest/$tag";
	    if (! -d $dest)
	    {
		mkdir($dest);
	    }
	}
	if (open($fh, ">", "$dest/$name"))
	{
	    $self->{file} = "$dest/$name";
	    $self->{dest} = $fh;
	}
	else
	{
	    warn "Cannot open log file $dest/$name: $!";
	}
    }
    else
    {
	#
	# Log to string.
	#
	my $stderr;
	$self->{dest} = \$stderr;
    }
    
    bless $self, $class;

    for my $e (sort { $a cmp $b } keys %ENV)
    {
	$self->log_cmd($e, $ENV{$e});
    }
    return $self;
}

sub redirect
{
    my($self) = @_;
    if ($self->{dest})
    {
	return("2>", $self->{dest});
    }
    else
    {
	return ();
    }
}

sub redirect_both
{
    my($self) = @_;
    if ($self->{dest})
    {
	return(">&", $self->{dest});
    }
    else
    {
	return ();
    }
}

sub timestamp
{
    my($self) = @_;
    my ($t, $us) = gettimeofday;
    $us = sprintf("%06d", $us);
    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
    return $ts;
}

sub log
{
    my($self, $str) = @_;
    my $d = $self->{dest};
    my $ts = $self->timestamp();
    if (ref($d) eq 'SCALAR')
    {
	$$d .= "[$ts] " . $str . "\n";
	return 1;
    }
    elsif ($d)
    {
	print $d "[$ts] " . $str . "\n";
	return 1;
    }
    return 0;
}

sub log_cmd
{
    my($self, @cmd) = @_;
    my $d = $self->{dest};
    my $str;
    my $ts = $self->timestamp();
    if (ref($cmd[0]))
    {
	$str = join(" ", @{$cmd[0]});
    }
    else
    {
	$str = join(" ", @cmd);
    }
    if (ref($d) eq 'SCALAR')
    {
	$$d .= "[$ts] " . $str . "\n";
    }
    elsif ($d)
    {
	print $d "[$ts] " . $str . "\n";
    }
	 
}

sub dest
{
    my($self) = @_;
    return $self->{dest};
}

sub text_value
{
    my($self) = @_;
    if (ref($self->{dest}) eq 'SCALAR')
    {
	my $r = $self->{dest};
	return $$r;
    }
    else
    {
	return $self->{file};
    }
}


1;
