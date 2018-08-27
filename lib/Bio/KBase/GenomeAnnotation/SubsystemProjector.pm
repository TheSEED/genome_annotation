package Bio::KBase::GenomeAnnotation::SubsystemProjector;

#
# This module wraps up the data required for the SEEDtk subsystem
# projector as well as the postprocessing required for use in the
# RASTtk genome annotation code.
#
# It stores lookup data in memory as loaded from flat files. It's mildly
# intended to be loaded at data curation time and exported via something
# like Sereal::Encode and reloaded from the file later when used.
#

use strict;
use SubsystemProjector;
use JSON::XS;
use File::Slurp;
use Encode;
use GenomeTypeObject;
use Data::Dumper;
use RoleParse;

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(projector role_file variant_file variant_code_map ref_ss));

sub new
{
    my($class, $role_file, $variant_file) = @_;

    my $proj = SubsystemProjector->new($role_file, $variant_file);
    my $self = {
	projector => $proj,
	role_file => $role_file,
	variant_file => $variant_file,
	variant_code_map => {},
	ref_ss => {},
    };
    return bless $self, $class;
}

#
# Project subsystems.
#
# Input hash may contain the following parameters:
#
# genome_object: project from the data in the given genome object
# genome_file: project from the data in the given file, assumed to be a genome object encoded as JSON
#
#
# Returns a subsystem list as defined in GenomeAnnotation.spec.
#

sub project_subsystems
{
    my($self, $input) = @_;

    ref($input) eq 'HASH' or die "Bio::KBase::GenomeAnnotation::SubsystemProjector::project_subsystems: parameter must be a hash reference";

    my %assigns;
    if (my $g = $input->{genome_object})
    {
	for my $f ($g->features())
	{
	    $assigns{$f->{id}} = $f->{function};
	}
    }
    elsif (my $f = $input->{genome_file})
    {
	my $g = GenomeTypeObject->new({file => $f});
	for my $f ($g->features())
	{
	    $assigns{$f->{id}} = $f->{function};
	}
    }
    else
    {
	die "No input source defined in Bio::KBase::GenomeAnnotation::SubsystemProjector::project_subsystems";
    }

    #
    # Event tagging for the RASTtk pipeline.
    #
    my @event_tag;
    if ($input->{event_id})
    {
	@event_tag = (event_id => $input->{event_id});
    }

    my $res = $self->projector->Project(\%assigns);

    #
    # Postprocess the raw projection.
    #
    # Map variant codes.
    # Determine classification
    # Find role list, and generate fid lists in role order.
    # Rename any "--_gjo" subsystems to lose that suffix.
    #

    my @mapped;
    while (my($ss_name, $ss) = each %$res)
    {
	my $ref_ss = $self->ref_ss->{$ss_name};
	my $ref_roles = $ref_ss->{role_name};
	my %ref_role_index = map { $ref_roles->[$_], $_ } 0..$#$ref_roles;
	my %ref_role_index_fuzzy = map {
	    my($text, $ec, $tc, $hypo) = RoleParse::Parse($ref_roles->[$_]);
	    (RoleParse::Normalize($text), $_);
	} 0..$#$ref_roles;

	my($vc, $roles) = @$ss;
	my $mapped_vc = $self->variant_code_map->{$vc};

	my $mapped_name = $ss_name;
	$mapped_name =~ s/_--gjo$//;

	my @rlist = map { { role_id => $_, features => [] } } @$ref_roles;

	for my $r (@$roles)
	{
	    my($role, $fid) = @$r;
	    my $ridx = $ref_role_index{$role};
	    if (!defined($ridx))
	    {
		# try a fuzzy match
		my($text, $ec, $tc, $hypo) = RoleParse::Parse($role);
		$ridx = $ref_role_index_fuzzy{RoleParse::Normalize($text)};
		if (defined($ridx))
		{
		    print STDERR "Fuzzy role match for '$role' => '$text' => '$ref_roles->[$ridx]' in $ss_name\n";
		}
		else
		{
		    die "Cannot find reference index for $role in $ss_name " . Dumper(\%ref_role_index);
		}
	    }
	    push(@{$rlist[$ridx]->{features}}, $fid);
	}
	next if $mapped_vc eq 'delete';
	push @mapped, {
	    name => $mapped_name,
	    classification => [ $ref_ss->{superclass}, $ref_ss->{class}, $ref_ss->{subclass} ],
	    variant_code => $mapped_vc,
	    role_bindings => \@rlist,
	    @event_tag,
	};
    }

    return \@mapped;
}

sub load_reference_subsystems
{
    my($self, $ref_ss_json) = @_;

    my $ref_ss = $self->ref_ss;
    my $json = JSON::XS->new;
    my $txt = scalar read_file($ref_ss_json);
    my $ref;
    eval {
	$ref = $json->decode($txt);
    };
    if ($@ =~ /malformed UTF-8/)
    {
	#
	# Try decoding as cp1252 since we have word copy and paste stuff in here.
	#
	
	eval {
	    $json->utf(0);
	    $ref = $json->decode(encode('UTF-8', decode('cp1252', $txt)));
	};
    }
    if ($@)
    {
	die "Cannot parse $ref_ss_json: $@";
    }
    for my $ent (@$ref)
    {
	$ref_ss->{$ent->{subsystem_id}} = $ent;
    }
}

sub load_variant_codes
{
    my($self, $variant_codes_file) = @_;

    open(my $fh, "<", $variant_codes_file) or die "Cannot read $variant_codes_file: $!";
    while (<$fh>)
    {
	chomp;
	my($vc, $mapped_vc) = split(/\t/);
	$self->variant_code_map->{$vc} = $mapped_vc;
    }
    close($fh);
}
1;
