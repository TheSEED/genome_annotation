#!/usr/bin/env perl

=head1 NAME

test-annotation-funcdef - Test a single GenomeAnnotation funcdef on a genome

=head1 SYNOPSIS

    test-annotation-funcdef --funcdef <funcdef_name> --genome <genome.json> [options]

    # Examples:
    test-annotation-funcdef --funcdef call_features_CDS_prodigal --genome genome.json
    test-annotation-funcdef --funcdef annotate_proteins_kmer_v2 --genome genome.json --output result.json
    test-annotation-funcdef --funcdef call_features_rRNA_SEED --genome genome.json --params '["SSU", "LSU"]'

=head1 DESCRIPTION

This script tests a single funcdef from the GenomeAnnotation.spec by:

1. Reading a JSON-formatted genome typed object
2. Instantiating the GenomeAnnotationImpl
3. Invoking the specified funcdef method
4. Returning the modified genome object

=head1 OPTIONS

=over 4

=item --funcdef <name>

The name of the funcdef to call (e.g., call_features_CDS_prodigal, annotate_proteins_kmer_v2)

=item --genome <file>

Path to a JSON file containing a genome typed object

=item --output <file>

Optional output file for the resulting genome (default: stdout)

=item --params <json>

Optional JSON-formatted parameters to pass to the funcdef (for methods that take extra params)

=item --list

List all available funcdefs that take a genomeTO and return a genomeTO

=item --help

Show this help message

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use JSON::XS;
use Data::Dumper;
use File::Slurp;

# Add lib paths for the genome_annotation module
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../../p3_core/lib";

use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use GenomeTypeObject;

# Funcdefs that transform genomeTO -> genomeTO
# Format: method_name => { params => param_type_name or undef, description => text }
my %FUNCDEFS = (
    # No extra params (just genomeTO)
    'annotate_genome'                   => { params => undef, desc => 'Run full annotation pipeline' },
    'assign_functions_to_CDSs'          => { params => undef, desc => 'Assign functions to CDS features' },
    'call_selenoproteins'               => { params => undef, desc => 'Call selenoprotein features' },
    'call_pyrrolysoproteins'            => { params => undef, desc => 'Call pyrrolysine protein features' },
    'call_features_selenoprotein'       => { params => undef, desc => 'Call selenoprotein features (alias)' },
    'call_features_pyrrolysoprotein'    => { params => undef, desc => 'Call pyrrolysine features (alias)' },
    'call_features_insertion_sequences' => { params => undef, desc => 'Call insertion sequence features' },
    'call_features_tRNA_trnascan'       => { params => undef, desc => 'Call tRNA features using tRNAscan-SE' },
    'call_RNAs'                         => { params => undef, desc => 'Call all RNA features' },
    'call_features_CDS_prodigal'        => { params => undef, desc => 'Call CDS features using Prodigal' },
    'call_features_CDS_genemark'        => { params => undef, desc => 'Call CDS features using GeneMarkS' },
    'call_features_CDS_phanotate'       => { params => undef, desc => 'Call CDS features using PHANOTATE (phage)' },
    'call_features_CDS_FragGeneScan'    => { params => undef, desc => 'Call CDS features using FragGeneScan' },
    'call_features_prophage_phispy'     => { params => undef, desc => 'Call prophage features using PhiSpy' },
    'call_features_strep_suis_repeat'   => { params => undef, desc => 'Call Streptococcus suis repeat features' },
    'call_features_strep_pneumo_repeat' => { params => undef, desc => 'Call Streptococcus pneumoniae repeat features' },
    'call_features_crispr'              => { params => undef, desc => 'Call CRISPR features' },
    'translate_untranslated_proteins'   => { params => undef, desc => 'Generate protein translations for untranslated CDS' },
    'annotate_special_proteins'         => { params => undef, desc => 'Annotate special proteins (AMR, virulence, etc.)' },
    'annotate_special_proteins_v2'      => { params => undef, desc => 'Annotate special proteins v2' },
    'annotate_families_figfam_v1'       => { params => undef, desc => 'Assign FIGfam protein families' },
    'annotate_families_patric'          => { params => undef, desc => 'Assign PATRIC protein families' },
    'annotate_families_patric_viral'    => { params => undef, desc => 'Assign PATRIC viral protein families' },
    'annotate_null_to_hypothetical'     => { params => undef, desc => 'Label unannotated proteins as hypothetical' },
    'remove_genbank_features'           => { params => undef, desc => 'Remove GenBank-imported features' },
    'annotate_strain_type_MLST'         => { params => undef, desc => 'Perform MLST strain typing' },
    'annotate_strain_type_MLST_v2'      => { params => undef, desc => 'Perform MLST strain typing v2' },
    'annotate_strain_type_cgMLST'       => { params => undef, desc => 'Perform cgMLST strain typing' },
    'annotate_proteins'                 => { params => undef, desc => 'Annotate all proteins' },
    'find_close_neighbors'              => { params => undef, desc => 'Find closely related genomes' },
    'renumber_features'                 => { params => undef, desc => 'Renumber feature IDs' },
    'classify_amr'                      => { params => undef, desc => 'Classify AMR genes' },
    'classify_amr_v2'                   => { params => undef, desc => 'Classify AMR genes v2' },
    'project_subsystems'                => { params => undef, desc => 'Project metabolic subsystems' },
    'compute_genome_quality_control'    => { params => undef, desc => 'Compute genome quality metrics' },
    'import_sra_metadata'               => { params => undef, desc => 'Import SRA metadata' },
    'compute_sars2_variation'           => { params => undef, desc => 'Compute SARS-CoV-2 variation' },

    # With parameters
    'call_features_rRNA_SEED'           => { params => 'list<rna_type>', desc => 'Call rRNA features (types: SSU, LSU, 5S)' },
    'call_features_CDS_glimmer3'        => { params => 'glimmer3_parameters', desc => 'Call CDS using Glimmer3' },
    'call_features_lowvan'              => { params => 'lowvan_parameters', desc => 'Call features using LoVaN' },
    'call_features_vigor4'              => { params => 'vigor4_parameters', desc => 'Call features using VIGOR4' },
    'call_features_vipr_mat_peptide'    => { params => 'vipr_mat_peptide_parameters', desc => 'Call VIPR mat peptide features' },
    'call_features_repeat_region_SEED'  => { params => 'repeat_region_SEED_parameters', desc => 'Call repeat region features' },
    'call_features_assembly_gap'        => { params => 'assembly_gap_parameters', desc => 'Call assembly gap features' },
    'call_features_CDS_SEED_projection' => { params => 'SEED_projection_parameters', desc => 'Call CDS using SEED projection' },
    'call_features_ProtoCDS_kmer_v1'    => { params => 'kmer_v1_parameters', desc => 'Call proto-CDS using kmer v1' },
    'call_features_ProtoCDS_kmer_v2'    => { params => 'kmer_v2_parameters', desc => 'Call proto-CDS using kmer v2' },
    'prune_invalid_CDS_features'        => { params => 'prune_invalid_CDS_features_parameters', desc => 'Prune invalid CDS features' },
    'annotate_proteins_similarity'      => { params => 'similarity_parameters', desc => 'Annotate proteins by BLAST similarity' },
    'annotate_proteins_phage'           => { params => 'phage_parameters', desc => 'Annotate proteins using phage DB' },
    'annotate_proteins_kmer_v1'         => { params => 'kmer_v1_parameters', desc => 'Annotate proteins using kmer v1' },
    'annotate_proteins_kmer_v2'         => { params => 'kmer_v2_parameters', desc => 'Annotate proteins using kmer v2' },
    'resolve_overlapping_features'      => { params => 'resolve_overlapping_features_parameters', desc => 'Resolve overlapping features' },
    'propagate_genbank_feature_metadata'=> { params => 'propagate_genbank_feature_metadata_parameters', desc => 'Propagate GenBank metadata' },
    'evaluate_genome'                   => { params => 'evaluate_genome_parameters', desc => 'Evaluate genome quality' },

    # Special cases with multiple extra params
    'call_features_scan_for_matches'    => { params => 'SPECIAL', desc => 'Scan for pattern matches (pattern, feature_type)' },
    'set_metadata'                      => { params => 'SPECIAL', desc => 'Set genome metadata' },
    'add_contigs'                       => { params => 'SPECIAL', desc => 'Add contigs to genome' },
    'add_features'                      => { params => 'SPECIAL', desc => 'Add features to genome' },
    'run_pipeline'                      => { params => 'SPECIAL', desc => 'Run annotation workflow/pipeline' },
);

my ($funcdef, $genome_file, $output_file, $params_json, $list_funcdefs, $help);

GetOptions(
    'funcdef=s'  => \$funcdef,
    'genome=s'   => \$genome_file,
    'output=s'   => \$output_file,
    'params=s'   => \$params_json,
    'list'       => \$list_funcdefs,
    'help'       => \$help,
) or pod2usage(2);

pod2usage(1) if $help;

if ($list_funcdefs) {
    print "Available funcdefs (genomeTO -> genomeTO):\n\n";
    print sprintf("  %-40s %s\n", "FUNCDEF", "DESCRIPTION");
    print "  " . ("-" x 80) . "\n";
    for my $name (sort keys %FUNCDEFS) {
        my $info = $FUNCDEFS{$name};
        my $param_info = $info->{params} ? " [requires params]" : "";
        print sprintf("  %-40s %s%s\n", $name, $info->{desc}, $param_info);
    }
    exit 0;
}

unless ($funcdef && $genome_file) {
    pod2usage("Error: --funcdef and --genome are required\n");
}

unless (exists $FUNCDEFS{$funcdef}) {
    die "Error: Unknown funcdef '$funcdef'. Use --list to see available funcdefs.\n";
}

if ($FUNCDEFS{$funcdef}{params} eq 'SPECIAL') {
    die "Error: funcdef '$funcdef' requires special handling not yet implemented.\n" .
        "       This funcdef takes multiple parameters beyond the genome object.\n";
}

# Read the genome JSON
print STDERR "Reading genome from $genome_file...\n";
my $genome_json = read_file($genome_file);
my $coder = JSON::XS->new->ascii->pretty->allow_nonref;
my $genome_raw = $coder->decode($genome_json);

# Initialize GenomeTypeObject with indexes
print STDERR "Initializing GenomeTypeObject...\n";
my $genome = GenomeTypeObject->initialize($genome_raw);

# Create the implementation object
print STDERR "Creating GenomeAnnotationImpl...\n";
my $impl = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();

# Set up fake context (required by some methods)
$Bio::KBase::GenomeAnnotation::Service::CallContext = {
    user_id => $ENV{USER} || 'test_user',
    authenticated => 1,
    token => $ENV{KB_AUTH_TOKEN} || '',
};

# Parse params if provided
my @extra_params;
if ($params_json) {
    my $params = $coder->decode($params_json);
    @extra_params = ref($params) eq 'ARRAY' ? @$params : ($params);
}
elsif ($FUNCDEFS{$funcdef}{params}) {
    # Provide empty params hash/array for methods that expect them
    my $param_type = $FUNCDEFS{$funcdef}{params};
    if ($param_type =~ /^list</) {
        @extra_params = ([]);
    } else {
        @extra_params = ({});
    }
    print STDERR "Note: Using empty default params for $funcdef\n";
}

# Invoke the funcdef
print STDERR "Calling $funcdef...\n";
my $result;
eval {
    $result = $impl->$funcdef($genome, @extra_params);
};

if ($@) {
    die "Error calling $funcdef: $@\n";
}

# Prepare result for output (strips internal indexes)
if (ref($result) && ref($result) ne 'HASH') {
    $result = $result->prepare_for_return();
}

# Output the result
my $output_json = $coder->encode($result);

if ($output_file) {
    write_file($output_file, $output_json);
    print STDERR "Result written to $output_file\n";
} else {
    print $output_json;
}

print STDERR "Done.\n";
