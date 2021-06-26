/*
  API Access to the Genome Annotation Service.

  Provides support for gene calling, functional annotation, re-annotation. Use to extract annotation in
formation about an existing genome, or to create new annotations.

 */
module GenomeAnnotation
{

    /*
     * This is a handle service handle object, used for by-reference
     * passing of data files.
     */
    typedef structure {
    string file_name;
    string id;
    string type;
    string url;
    string remote_md5;
    string remote_sha1;
    } Handle;


    typedef int bool;
    typedef string md5;
    typedef list<md5> md5s;
    typedef string genome_id;
    typedef string feature_id;
    typedef string contig_id;
    typedef string feature_type;

    /* A region of DNA is maintained as a tuple of four components:

        the contig
        the beginning position (from 1)
        the strand
        the length

       We often speak of "a region".  By "location", we mean a sequence
       of regions from the same genome (perhaps from distinct contigs).

       Strand is either '+' or '-'.
        */
    typedef tuple<contig_id, int begin, string strand,int length> region_of_dna;

    /*
    a "location" refers to a sequence of regions
    */
    typedef list<region_of_dna> location;

    typedef string analysis_event_id;
    typedef structure {
    analysis_event_id id;
    string tool_name;
    float execution_time;
    list<string> parameters;
    string hostname;
    } analysis_event;

    typedef structure {
    string job_id;
    string start_time;
    string completion_time;
    float elapsed_time;
    string app_name;
    mapping<string, string> parameters;
    mapping<string, string> attributes;
    } job_statistics;

    typedef tuple<string comment, string annotator, float annotation_time, analysis_event_id> annotation;

    typedef structure {
    bool truncated_begin;
    bool truncated_end;
    /* Is this a real feature? */
    float existence_confidence;

    bool frameshifted;
    bool selenoprotein;
    bool pyrrolysylprotein;

    /*
     * List of rules that govern the overlap removal procedure for
     * this feature. We don't yet have a strict definition for this but
     * the notion is that this will consiste of entries of the form
     * +feature-type which will allow overlap with the given feature type;
     * -feature-type which will disallow overlap with the given feature type.
     */
    list<string> overlap_rules;

    /*
     * The numeric priority of this feature's right to exist. Specialty
     * tools will give the features they create a high priority; more generic
     * tools will give their features a lower priority. The overlap removal procedure
     * will use this priority to determine which of a set of overlapping features
     * should be removed.
     *
     * The intent is that a change of 1 in the priority value represents a factor of 2 in
     * preference.
     */
    float existence_priority;

    float hit_count;
    float weighted_hit_count;
    float genemark_score;
    } feature_quality_measure;

    /*
     * A protein family assignment notes the assignment of the given feature
     * to a protein family. db is the name of the protein family database
     * (e.g. FIGfam, GPF for GlobalPatricFam, LPF for LocalPatricFam, etc.)
     */

    typedef tuple <string db, string id, string function, string db_version> protein_family_assignment;

    /*
     * A similarity association notes the BLAST-computed association
     * between this feature and a given protein database.
     */

    typedef tuple <string source, string source_id,
    float query_coverage, float subject_coverage, float identity, float e_value>
    similarity_association;

    /* A proposed function records an assertion of the function of a feature.
     * A feature may have multiple proposed functions. A tool downstream of the
     * tools that propose functions may determine based on the asserted proposals
     * which function should be the assigned function for the feature.
     */
    typedef structure {
    string id;
    string function;
    string user;
    float score;
    analysis_event_id event_id;
    int timestamp;
    } proposed_function;

    typedef structure {
        string genbank_type;
        string genbank_location;
        mapping<string qualifier, list<string>> values;
    } genbank_feature;

    typedef structure {
    list<string> accession;
    list<string> comment;
    string date;
    list<string> dblink;
    list<string> dbsource;
    string definition;
    string division;
    string geometry;
    int gi;
    list<string> keywords;
    string locus;
    string organism;
    string origin;
    list<mapping<string, string>> references;
    string source;
    list<string> taxonomy;
    list<string> version;
    } genbank_locus;

    typedef tuple <string, int, double> feature_coupling;

    /* A feature object represents a feature on the genome. It contains
       the location on the contig with a type, the translation if it
       represents a protein, associated aliases, etc. It also contains
       information gathered during the annotation process that is involved
       in stages that perform overlap removal, quality testing, etc.
    */
    typedef structure {
    feature_id id;
    location location;
    feature_type type;
    string function;
    /*
     * The function_id refers to the particular proposed function that was chosen
     * for this feature.
     */
    string function_id;
    string protein_translation;
    list<string> aliases;
    list<tuple<string source, string alias>> alias_pairs;
    list<annotation> annotations;
    feature_quality_measure quality;
    analysis_event_id feature_creation_event;
    list<protein_family_assignment> family_assignments;
    list<similarity_association> similarity_associations;
    list<proposed_function> proposed_functions;

    string genbank_type;
    genbank_feature genbank_feature;

    list<tuple<string id, string description>> ec_numbers;
    list<tuple<string id, string description>> go_terms;
    list<tuple<string id, string description>> pathways;
    list<feature_coupling> couplings;
    } feature;

    /* Data for DNA contig */
    typedef structure {
        contig_id id;
        string dna;
        int genetic_code;
        string cell_compartment;
        string replicon_type;
        /* circular / linear */
        string replicon_geometry;
        bool complete;
        genbank_locus genbank_locus;
        string original_id;
    } contig;

    typedef structure {
    genome_id genome;
    string genome_name;
    float closeness_measure;
    string analysis_method;
    } close_genome;

    typedef structure
    {
    float frameshift_error_rate;
    float sequence_error_rate;
    structure {
        mapping<string, string> role_map;
        mapping<string, list<string>> role_fids;
        mapping<string, tuple<int predicted, int actual>> consistency_roles;
        mapping<string, tuple<int predicted, int actual>> completeness_roles;
        int predicted_roles;
        int under_present;
        int over_present;
        int consistency_checked;
        int completeness_checked;
    } problematic_roles_report;
    string eval_version;
    bool hasSsuRRna;
    string dna_md5;
    float coarse_consistency;
    float fine_consistency;
    float completeness;
    float contamination;
    string completeness_group;
    structure {
        int N50;
        int N70;
        int N90;
        int L50;
        int L70;
        int L90;
        int totlen;
        int complete;
    } genome_metrics;

    int genome_length;
    float gc_content;
    int chromosomes;
    int plasmids;
    int contigs;

    int contig_ambig_count;
    float contig_ambig_fraction;
    int contig_longest_ambig_run;

    string genome_status;

    structure {
        int cds;
        int partial_cds;
        int rRNA;
        int tRNA;
        int misc_RNA;
        int repeat_region;
    } feature_summary;

    structure {
        int hypothetical;
        int function_assignment;
        int plfam_assignment;
        int pgfam_assignment;
        int ec_assignment;
        int go_assignment;
        int pathway_assignment;
    } protein_summary;

    mapping<string, int> specialty_gene_summary;

    list<tuple<feature_id id, string gene_name, string function, string amr_classification>> amr_genes;
    list<tuple<string amr_classification, list<string> gene_names>> amr_gene_summary;

    mapping<string superclass, structure { int genes; int subsystems; } > subsystem_summary;

    float cds_ratio;
    float hypothetical_cds_ratio;
    float partial_cds_ratio;
    float plfam_cds_ratio;
    float pgfam_cds_ratio;

    list<string> genome_quality_flags;
    string genome_quality;

    } genome_quality_measure;

    typedef structure
    {
    string typing_method;
    string database;
    string tag;
    analysis_event_id event_id;
    } strain_type;

    typedef structure
    {
    string name;
    string version;
    string description;
    string comment;
    list<string> antibiotics;
    float accuracy;
    float area_under_roc_curve;
    float f1_score;
    string sources;
    float cumulative_adaboost_value;
    string sensitivity;
    analysis_event_id event_id;
    list<tuple<feature_id id, float alpha, int round, string function>> features;
    } classifier;

    typedef structure {
    string role_id;
    list<feature_id> features;
    } role_binding;

    typedef structure {
    string name;
    tuple<string superclass, string class, string subclass> classification;
    string variant_code;
    list<role_binding> role_bindings;
    analysis_event_id event_id;
    } subsystem_data;

    /* Variant support */

    typedef structure {
        int pos;
        string ref;
        string alt;
        float freq;
        string feature_pos;
        string ref_aa;
        string ref_codon;
        string alt_aa;
        string alt_codon;
    } snp;
    typedef structure {
        string reference;
        string gene;
        int frame;
        list<snp> snps;
    } variant;
    typedef structure {
        string tool;
        mapping<string, string> tool_metadata;
        list<variant> variants;
        string lineage;
        float probability;
        string status;
        string notes;
    } computed_variant;



    /* All of the information about particular genome */
    typedef structure {
    genome_id id;
    string scientific_name;
    string domain;
    int genetic_code;
    string source;
    string source_id;
    string taxonomy;
    int ncbi_taxonomy_id;
    list<tuple<string taxon_name, int taxon_id, string taxon_rank>> ncbi_lineage;
    string ncbi_genus;
    string ncbi_species;
    string owner;
    string home;

    genome_quality_measure quality;

    list<contig> contigs;
    Handle contigs_handle;

    list<feature> features;

    list<close_genome> close_genomes;

    list <analysis_event> analysis_events;

    list<strain_type> typing;
    list<classifier> classifications;

    list<subsystem_data> subsystems;

    structure {
        job_statistics assembly;
        job_statistics annotation;
    } job_data;

    list<mapping<string key, string value>> sra_metadata;

    list<computed_variant> computed_variants;




    } genomeTO;


    /*
     * Genome metadata. We use this structure to define common metadata
     * settings used in the API calls below. It is possible this data should
     * have been separated in this way in the genome object itself, but there
     * is an extant body of code that assumes the current structure of the genome
     * object.
     */

    typedef structure
    {
    genome_id id;
    string scientific_name;
    string domain;
    int genetic_code;
    string source;
    string source_id;
    int ncbi_taxonomy_id;
    string taxonomy;
    string owner;
    } genome_metadata;

    typedef string subsystem;
    typedef string variant;
    typedef tuple<subsystem,variant> variant_of_subsystem;
    typedef list<variant_of_subsystem> variant_subsystem_pairs;
    typedef string fid;
    typedef string role;
    typedef string function;
    typedef tuple<fid,role> fid_role_pair;
    typedef list<fid_role_pair> fid_role_pairs;
    typedef tuple<fid,function> fid_function_pair;
    typedef list<fid_function_pair> fid_function_pairs;

    /* Metabolic reconstruction
       represents the set of subsystems that we infer are present in this genome
    */
    typedef structure {
    variant_subsystem_pairs subsystems;
    fid_role_pairs bindings;
    fid_function_pairs assignments;
    } reconstructionTO;

    typedef tuple<fid,md5,location,function> fid_data_tuple;
    typedef list<fid_data_tuple> fid_data_tuples;

    /*
     * Given one or more Central Store genome IDs, convert them into genome objects.
     */
    funcdef genome_ids_to_genomes(list<genome_id> ids) returns (list<genomeTO> genomes);

    /*
     * Create a new genome object and assign metadata.
     */
    funcdef create_genome(genome_metadata metadata) returns (genomeTO genome);

    /*
     * Create a new genome object from one or more genbank files.
     */
    funcdef create_genome_from_genbank(string gb_data) returns (genomeTO genome);

    /*
     * Create a new genome object based on data from the SEED project.
     */
    funcdef create_genome_from_SEED(string genome_id) returns (genomeTO genome);

    /*
     * Create a new genome object based on a RAST genome.
     */
    funcdef create_genome_from_RAST(string genome_or_job_id) returns (genomeTO genome); /*  authentication optional; */

    /*
     * Modify genome metadata.
     */
    funcdef set_metadata(genomeTO genome_in, genome_metadata metadata) returns (genomeTO genome_out);

    /*
     * Add a set of contigs to the genome object.
     */
    funcdef add_contigs(genomeTO genome_in, list<contig> contigs) returns (genomeTO genome_out);

    /*
     * Add a set of contigs to the genome object, loading the contigs
     * from the given handle service handle.
     */
    funcdef add_contigs_from_handle(genomeTO genome_in, list<contig> contigs) returns (genomeTO genome_out);

    /*
     * Import SRA metadata from initial assembly, if present.
     */
    funcdef import_sra_metadata(genomeTO genome_in) returns (genomeTO genome_out);

    /*
     * Compute SARS2 variation data.
     */
    funcdef compute_sars2_variation(genomeTO genome_in) returns ( genomeTO genome_out);

    /*
     * This tuple defines a compact form for defining features to be batch-loaded
     * into a genome object.
     */
    typedef tuple <string id, string location, string feature_type, string function, string aliases> compact_feature;

    /*
     * Add a set of features in tabular form.
     */
    funcdef add_features(genomeTO genome_in, list<compact_feature> features) returns (genomeTO genome_out);

    funcdef genomeTO_to_reconstructionTO (genomeTO) returns (reconstructionTO);
    funcdef genomeTO_to_feature_data (genomeTO) returns (fid_data_tuples);
    funcdef reconstructionTO_to_roles (reconstructionTO) returns (list<role>);
    funcdef reconstructionTO_to_subsystems(reconstructionTO) returns (variant_subsystem_pairs);

    /*
     * Given a genome object populated with contig data, perform gene calling
     * and functional annotation and return the annotated genome.
     */
    funcdef assign_functions_to_CDSs(genomeTO) returns (genomeTO);
    funcdef annotate_genome(genomeTO) returns (genomeTO);

    funcdef call_selenoproteins(genomeTO) returns (genomeTO);
    funcdef call_pyrrolysoproteins(genomeTO) returns (genomeTO);

    /*
     * Given a genome typed object, call selenoprotein features.
     */
    funcdef call_features_selenoprotein(genomeTO) returns (genomeTO);

    /*
     * Given a genome typed object, call pyrrolysoprotein features.
     */
    funcdef call_features_pyrrolysoprotein(genomeTO) returns (genomeTO);

    /*
     * Given a genome typed object, call insertion sequences.
     */
    funcdef call_features_insertion_sequences(genomeTO) returns (genomeTO);

    /* [ validate.enum("5S", "SSU", "LSU", "ALL") ] */
    typedef string rna_type;

    /*
     * Given a genome typed object, find instances of ribosomal RNAs in
     * the genome.
     *
     * The types parameter is used to select the types of RNAs to
     * call. It is a list of strings where each value is one of
     *
     *    "5S"
     *    "SSU"
     *    "LSU"
     *
     * or "ALL" to choose all available rRNA types.
     */
    funcdef call_features_rRNA_SEED(genomeTO genome_in, list<rna_type> types) returns (genomeTO genome_out);

    /*
     * Given a genome typed object, find instances of tRNAs in
     * the genome.
     */
    funcdef call_features_tRNA_trnascan(genomeTO genome_in) returns (genomeTO genome_out);

    /*
     * Given a genome typed object, find instances of all RNAs we currently
     * have support for detecting.
     */
    funcdef call_RNAs(genomeTO genome_in) returns (genomeTO genome_out);

    typedef structure
    {
    string reference_name;
    bool remove_existing_features;
    } vigor4_parameters;

    funcdef call_features_vigor4(genomeTO, vigor4_parameters params) returns (genomeTO);

    typedef structure
    {
    bool remove_existing_features;
    } vipr_mat_peptide_parameters;

    funcdef call_features_vipr_mat_peptide(genomeTO, vipr_mat_peptide_parameters params) returns (genomeTO);

    typedef structure
    {
    int min_training_len;
    } glimmer3_parameters;

    funcdef call_features_CDS_glimmer3(genomeTO, glimmer3_parameters params) returns (genomeTO);

    funcdef call_features_CDS_prodigal(genomeTO) returns (genomeTO);
    funcdef call_features_CDS_genemark(genomeTO) returns (genomeTO);
    funcdef call_features_CDS_phanotate(genomeTO) returns (genomeTO);

    typedef structure
    {
    int minimum_contig_length;
    float max_homopolymer_frequency;
    float max_dna_in_translation;
    } prune_invalid_CDS_features_parameters;
    funcdef prune_invalid_CDS_features(genomeTO genome_in, prune_invalid_CDS_features_parameters params) returns (genomeTO genome_out);

    typedef structure
    {
    string reference_database;
    string reference_id;
    int kmer_size;
    } SEED_projection_parameters;

    funcdef call_features_CDS_SEED_projection(genomeTO, SEED_projection_parameters params) returns (genomeTO);
    funcdef call_features_CDS_FragGeneScan(genomeTO) returns (genomeTO);

    typedef structure
    {
    float min_identity;
    int min_length;
    } repeat_region_SEED_parameters;
    funcdef call_features_repeat_region_SEED(genomeTO genome_in, repeat_region_SEED_parameters params) returns (genomeTO genome_out);

    funcdef call_features_prophage_phispy(genomeTO genome_in) returns (genomeTO genome_out);

    funcdef call_features_scan_for_matches(genomeTO genome_in, string pattern, string feature_type) returns (genomeTO genome_out);

    typedef structure
    {
    int min_gap_length;
    int monopolymer_repeat_length;
    } assembly_gap_parameters;

    /*
    * Given a genome typed object, call gap features.
    * Gaps are known regions in the contig where the nucleotide sequence is not known
    * but where there is evidence that a run of DNA does exist joining the sequenced
    * data on either side of the gap.
    *
    * Gaps are currently called using one of two methods. Genomes that originated as
    * genbank files may have a CONTIGS entry that defines the contig and gap regions.
    * Genomes that do not have a CONTIGS entry are scanned for runs of "n" characters.
    */
    funcdef call_features_assembly_gap(genomeTO genome_in, assembly_gap_parameters params) returns (genomeTO genome_out);

    typedef structure
    {
    int tmp;
    } split_gap_spanning_features_params;

    funcdef split_gap_spanning_features(genomeTO genome_in, split_gap_spanning_features_params)
            returns (genomeTO genome_out);

    funcdef translate_untranslated_proteins(genomeTO genome_in) returns (genomeTO genome_out);

    typedef structure
    {
    int annotate_hypothetical_only;
    int annotate_null_only;
    } similarity_parameters;

    /*
     * Annotate based on similarity to annotation databases.
     */
    funcdef annotate_proteins_similarity(genomeTO, similarity_parameters params) returns (genomeTO);

    typedef structure
    {
    int annotate_hypothetical_only;
    int annotate_null_only;
    } phage_parameters;
    /*
     * Annotate based on similarity to the phage annotation daatabase.
     */
    funcdef annotate_proteins_phage(genomeTO, phage_parameters params) returns (genomeTO);

    typedef structure
    {
    int kmer_size;
    string dataset_name;
    int return_scores_for_all_proteins;
    int score_threshold;
    int hit_threshold;
    int sequential_hit_threshold;
    int detailed;
    int min_hits;
    int min_size;
    int max_gap;
    int annotate_hypothetical_only;
    int annotate_null_only;
    } kmer_v1_parameters;

    funcdef annotate_proteins_kmer_v1(genomeTO, kmer_v1_parameters params) returns (genomeTO);

    typedef structure {
    int min_hits;
    int max_gap;
    int annotate_hypothetical_only;
    int annotate_null_only;
    string kmer_data_directory;
    string kmer_service_url;
    } kmer_v2_parameters;

    funcdef annotate_proteins_kmer_v2(genomeTO genome_in, kmer_v2_parameters params) returns (genomeTO genome_out);

    typedef structure {
    int placeholder;
    } resolve_overlapping_features_parameters;

    funcdef resolve_overlapping_features(genomeTO genome_in, resolve_overlapping_features_parameters params) returns (genomeTO genome_out);

    typedef structure {
    float min_rna_pct_coverage;
    } propagate_genbank_feature_metadata_parameters;

    funcdef propagate_genbank_feature_metadata(genomeTO genome_in, propagate_genbank_feature_metadata_parameters params) returns (genomeTO genome_out);

    funcdef call_features_ProtoCDS_kmer_v1(genomeTO, kmer_v1_parameters params) returns (genomeTO);
    funcdef call_features_ProtoCDS_kmer_v2(genomeTO genome_in, kmer_v2_parameters params) returns (genomeTO genome_out);

    funcdef enumerate_special_protein_databases() returns (list<string> database_names);

    typedef tuple <
    string protein_id,
    string database_name,
    string database_id,
    string protein_coverage,
    string database_coverage,
    float identity,
    float p_value > special_protein_hit;
    funcdef compute_special_proteins(genomeTO genome_in, list<string> database_names) returns (list<special_protein_hit> results);

    funcdef annotate_special_proteins(genomeTO genome_in) returns (genomeTO genome_out);
    funcdef annotate_families_figfam_v1(genomeTO genome_in) returns (genomeTO genome_out);
    funcdef annotate_families_patric(genomeTO genome_in) returns (genomeTO genome_out);
    funcdef annotate_families_patric_viral(genomeTO genome_in) returns (genomeTO genome_out);
    funcdef annotate_null_to_hypothetical(genomeTO genome_in) returns (genomeTO genome_out);

    funcdef remove_genbank_features(genomeTO genome_in) returns (genomeTO genome_out);

    funcdef annotate_strain_type_MLST(genomeTO genome_in) returns (genomeTO genome_out);

    typedef tuple <
    string protein_id,
    string domain_id,
    float identity,
    int alignment_len,
    int mismatches,
    int gap_openings,
    int protein_start,
    int protein_end,
    int domain_start,
    int domain_end,
    float e_value,
    float bit_score,
    string accession,
    string short_name,
    string description,
    int pssm_length > cdd_hit;
    funcdef compute_cdd(genomeTO genome_in) returns (list<cdd_hit>);

    funcdef annotate_proteins(genomeTO) returns (genomeTO);

    /* Determine close genomes. */
    funcdef estimate_crude_phylogenetic_position_kmer(genomeTO) returns (string position_estimate);

    funcdef find_close_neighbors(genomeTO) returns (genomeTO);

    /*
     * Interface to Strep repeats and "boxes" tools
     */
    funcdef call_features_strep_suis_repeat(genomeTO) returns (genomeTO);
    funcdef call_features_strep_pneumo_repeat(genomeTO) returns (genomeTO);
    funcdef call_features_crispr(genomeTO genome_in) returns (genomeTO genome_out);

    funcdef update_functions(genomeTO genome_in, list<tuple<feature_id, string function>> functions, analysis_event event)
    returns (genomeTO genome_out);

    /*
     * Renumber features such that their identifiers are contiguous along contigs.
     *
     */
    funcdef renumber_features(genomeTO genome_in) returns (genomeTO genome_out);

    /*
     * Perform AMR classification.
     */
    funcdef classify_amr(genomeTO) returns (genomeTO);

    typedef structure {
    string reference_genome_id;
    } evaluate_genome_parameters;
    /*
     * Perform genome evaluation.
     */
    funcdef evaluate_genome(genomeTO genome_in, evaluate_genome_parameters params) returns (genomeTO genome_out);

    /*
     * Export genome typed object to one of the supported output formats:
     * genbank, embl, or gff.
     * If feature_types is a non-empty list, limit the output to the given
     * feature types.
     */
    funcdef export_genome(genomeTO genome_in, string format, list<string> feature_types) returns (string exported_data);

    /*
     * Enumerate the available classifiers. Returns the list of identifiers for
     * the classifiers.
     */
    funcdef enumerate_classifiers() returns (list<string>);

    /*
     * Query the groups included in the given classifier. This is a
     * mapping from the group name to the list of genome IDs included
     * in the group. Note that these are genome IDs native to the
     * system that created the classifier; currently these are
     * SEED genome IDs that may be translated using the
     * source IDs on the Genome entity.
     */
    funcdef query_classifier_groups(string classifier) returns(mapping<string group_id, list<genome_id>>);

    /*
     * Query the taxonomy strings that this classifier maps.
     */

    funcdef query_classifier_taxonomies(string classifier) returns(mapping<string group_id, string taxonomy>);

    /*
     * Classify a dataset, returning only the binned output.
     */
    funcdef classify_into_bins(string classifier, list<tuple<string id, string dna_data>> dna_input)
    returns(mapping<string group_id, int count>);

    /*
     * Classify a dataset, returning the binned output along with the raw assignments and the list of
     * sequences that were not assigned.
     */
    funcdef classify_full(string classifier, list<tuple<string id, string dna_data>> dna_input)
    returns(mapping<string group_id, int count>, string raw_output, list<string> unassigned);


    /*
     * Project subsystems.
     */
    funcdef project_subsystems(genomeTO genome_in) returns (genomeTO genome_out);

    typedef structure {
    string name;
    string condition;
    int failure_is_not_fatal;
    repeat_region_SEED_parameters repeat_region_SEED_parameters;
    glimmer3_parameters glimmer3_parameters;
    kmer_v1_parameters kmer_v1_parameters;
    kmer_v2_parameters kmer_v2_parameters;
    similarity_parameters similarity_parameters;
    } pipeline_stage;

    /*
     * Compute genome Quality Control scoring.
     */
    funcdef compute_genome_quality_control(genomeTO genome_in) returns (genomeTO genome_out);

    typedef structure
    {
    list<pipeline_stage> stages;
    } workflow;

    typedef structure
    {
    string id;
    string name;
    string description;
    workflow workflow;
    } recipe;

    funcdef default_workflow() returns (workflow);

    /*
     * Enumerate the loaded workflows. We always have a workflow named "default"; a
     * particular deployment of the genome annotation service may include additional workflows.
     */
    funcdef enumerate_recipes() returns (list<recipe> recipes);

    /*
     * Look up and return a particular named workflow.
     */
    funcdef find_recipe(string id) returns (recipe);

    /*
     * Return a workflow that includes all available stages. Not meant
     * (necessarily) for actual execution, but as a comprehensive list
     * of parts for users to use in assembling their own workflows.
     */
    funcdef complete_workflow_template() returns (workflow);
    funcdef run_pipeline(genomeTO genome_in, workflow workflow) returns (genomeTO genome_out);

    typedef structure
    {
    string genome_id;
    Handle data;
    string filename;
    } pipeline_batch_input;

    typedef structure
    {
    string genome_id;
    string status;
    string creation_date;
    string start_date;
    string completion_date;

    Handle stdout;
    Handle stderr;
    Handle output;

    string filename;
    } pipeline_batch_status_entry;

    typedef structure
    {
    string status;
    string submit_date;
    string start_date;
    string completion_date;

    list<pipeline_batch_status_entry> details;
    } pipeline_batch_status;

    funcdef pipeline_batch_start(list<pipeline_batch_input> genomes, workflow workflow)
    returns (string batch_id) authentication required;
    funcdef pipeline_batch_status(string batch_id)
    returns (pipeline_batch_status status) authentication required;
    funcdef pipeline_batch_enumerate_batches()
    returns (list<tuple<string batch_id, string submit_time>> batches) authentication required;
};
