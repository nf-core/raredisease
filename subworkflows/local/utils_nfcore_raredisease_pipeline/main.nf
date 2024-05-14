//
// Subworkflow with functionality specific to the nf-core/raredisease pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    Channel.fromSamplesheet("input")
        .tap { ch_original_input }
        .map { meta, fastq1, fastq2 -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2 ->
            new_meta = meta + [num_lanes:counts[meta.id],
                        read_group:"\'@RG\\tID:"+ fastq1.toString().split('/')[-1] + "\\tPL:" + params.platform.toUpperCase() + "\\tSM:" + meta.id + "\'"]
            if (!fastq2) {
                return [ new_meta + [ single_end:true ], [ fastq1 ] ]
            } else {
                return [ new_meta + [ single_end:false ], [ fastq1, fastq2 ] ]
            }
        }
        .tap{ ch_input_counts }
        .map { meta, fastqs -> fastqs }
        .reduce([:]) { counts, fastqs -> //get line number for each row to construct unique sample ids
            counts[fastqs] = counts.size() + 1
            return counts
        }
        .combine( ch_input_counts )
        .map { lineno, meta, fastqs -> //append line number to sampleid
            new_meta = meta + [id:meta.id+"_T"+lineno[fastqs]]
            return [ new_meta, fastqs ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def align_text           = []
    def variant_call_text    = []
    def repeat_call_text     = []
    def snv_annotation_text  = []
    def sv_annotation_text   = []
    def mt_annotation_text   = []
    def qc_bam_text          = []
    def me_calls_text        = []
    def me_annotation_text   = []
    def preprocessing_text   = []
    def other_citation_text  = []

    align_text = [
        params.aligner.equals("bwa") ? "BWA (Li, 2013)," :"",
        params.aligner.equals("bwamem2") ? "BWA-MEM2 (Vasimuddin et al., 2019)," : "",
        params.aligner.equals("sentieon") ? "Sentieon DNASeq (Kendig et al., 2019)," : "",
        params.aligner.equals("sentieon") ? "Sentieon Tools (Freed et al., 2017)," : ""
    ]
    variant_call_text = [
        params.variant_caller.equals("deepvariant") ? "DeepVariant (Poplin et al., 2018)," : "",
        params.variant_caller.equals("sentieon") ? "Sentieon DNAscope (Freed et al., 2022)," : "",
        params.skip_haplocheck ? "" : "Haplocheck (Weissensteiner et al., 2021),",
        "CNVnator (Abyzov et al., 2011),",
        "TIDDIT (Eisfeldt et al., 2017),",
        "Manta (Chen et al., 2016),",
        "GLnexus (Yun et al., 2021),",
        params.skip_eklipse ? "" : "eKLIPse (Goudenge et al., 2019),",
    ]
    repeat_call_text = [
        "ExpansionHunter (Dolzhenko et al., 2019),",
        "stranger (Nilsson & Magnusson, 2021),"
    ]
    if (!params.skip_snv_annotation) {
        snv_annotation_text = [
            "CADD (Rentzsch et al., 2019, 2021),",
            "Vcfanno (Pedersen et al., 2016),",
            "VEP (McLaren et al., 2016),",
            "Genmod (Magnusson et al., 2018),",
        ]
    }
    if (!params.skip_sv_annotation) {
        sv_annotation_text = [
            "SVDB (Eisfeldt et al., 2017),",
            "VEP (McLaren et al., 2016),",
            "Genmod (Magnusson et al., 2018),",
        ]
    }
    if (!params.skip_mt_annotation) {
        mt_annotation_text = [
            "CADD (Rentzsch et al., 2019, 2021),",
            "VEP (McLaren et al., 2016),",
            "Vcfanno (Pedersen et al., 2016),",
            "Hmtnote (Preste et al., 2019),",
            "HaploGrep2 (Weissensteiner et al., 2016),",
            "Genmod (Magnusson et al., 2018),",
        ]
    }
    if (!params.skip_me_annotation) {
        me_annotation_text = [
            "VEP (McLaren et al., 2016),",
            "SVDB (Eisfeldt et al., 2017),",
        ]
    }
    qc_bam_text = [
        "Picard (Broad Institute, 2023)",
        "Qualimap (Okonechnikov et al., 2016),",
        "TIDDIT (Eisfeldt et al., 2017),",
        "UCSC Bigwig and Bigbed (Kent et al., 2010),",
        "Mosdepth (Pedersen & Quinlan, 2018),",
    ]
    me_calls_text = [
        "SVDB (Eisfeldt et al., 2017),",
        "RetroSeq (Keane et al., 2013),",
    ]
    preprocessing_text = [
        params.skip_fastqc ? "" : "FastQC (Andrews 2010),",
        params.skip_fastp  ? "" : "Fastp (Chen, 2023),",
    ]
    other_citation_text = [
        "BCFtools (Danecek et al., 2021),",
        "GATK (McKenna et al., 2010),",
        "MultiQC (Ewels et al. 2016),",
        params.skip_peddy ? "" : "Peddy (Pedersen & Quinlan, 2017),",
        params.run_rtgvcfeval ? "RTG Tools (Cleary et al., 2015)," : "",
        "SAMtools (Li et al., 2009),",
        "SMNCopyNumberCaller (Chen et al., 2020),",
        "Tabix (Li, 2011)",
        "."
    ]

    def concat_text = align_text +
                        variant_call_text    +
                        repeat_call_text     +
                        snv_annotation_text  +
                        sv_annotation_text   +
                        mt_annotation_text   +
                        qc_bam_text          +
                        me_calls_text        +
                        me_annotation_text   +
                        preprocessing_text   +
                        other_citation_text

    def citation_text = [ "Tools used in the workflow included:" ] + concat_text.unique(false) { a, b -> a <=> b } - ""
    return citation_text.join(' ').trim()
}

def toolBibliographyText() {

    def align_text           = []
    def variant_call_text    = []
    def repeat_call_text     = []
    def snv_annotation_text  = []
    def sv_annotation_text   = []
    def mt_annotation_text   = []
    def qc_bam_text          = []
    def me_calls_text        = []
    def me_annotation_text   = []
    def preprocessing_text   = []
    def other_citation_text  = []

    align_text = [
        params.aligner.equals("bwa") ? "<li>Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM (arXiv:1303.3997). arXiv. http://arxiv.org/abs/1303.3997</li>" :"",
        params.aligner.equals("bwamem2") ? "<li>Vasimuddin, Md., Misra, S., Li, H., & Aluru, S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), 314–324. https://doi.org/10.1109/IPDPS.2019.00041</li>" : "",
        params.aligner.equals("sentieon") ? "<li>Kendig, K. I., Baheti, S., Bockol, M. A., Drucker, T. M., Hart, S. N., Heldenbrand, J. R., Hernaez, M., Hudson, M. E., Kalmbach, M. T., Klee, E. W., Mattson, N. R., Ross, C. A., Taschuk, M., Wieben, E. D., Wiepert, M., Wildman, D. E., & Mainzer, L. S. (2019). Sentieon DNASeq Variant Calling Workflow Demonstrates Strong Computational Performance and Accuracy. Frontiers in Genetics, 10, 736. https://doi.org/10.3389/fgene.2019.00736</li>" : "",
        params.aligner.equals("sentieon") ? "<li>Freed, D., Aldana, R., Weber, J. A., & Edwards, J. S. (2017). The Sentieon Genomics Tools—A fast and accurate solution to variant calling from next-generation sequence data (p. 115717). bioRxiv. https://doi.org/10.1101/115717</li>" : ""
    ]
    variant_call_text = [
        params.variant_caller.equals("deepvariant") ? "<li>Poplin, R., Chang, P.-C., Alexander, D., Schwartz, S., Colthurst, T., Ku, A., Newburger, D., Dijamco, J., Nguyen, N., Afshar, P. T., Gross, S. S., Dorfman, L., McLean, C. Y., & DePristo, M. A. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983–987. https://doi.org/10.1038/nbt.4235</li>" : "",
        params.variant_caller.equals("sentieon") ? "<li>Freed, D., Pan, R., Chen, H., Li, Z., Hu, J., & Aldana, R. (2022). DNAscope: High accuracy small variant calling using machine learning [Preprint]. Bioinformatics. https://doi.org/10.1101/2022.05.20.492556</li>" : "",
        params.skip_haplocheck ? "" : "<li>Weissensteiner, H., Forer, L., Fendt, L., Kheirkhah, A., Salas, A., Kronenberg, F., & Schoenherr, S. (2021). Contamination detection in sequencing studies using the mitochondrial phylogeny. Genome Research, 31(2), 309–316. https://doi.org/10.1101/gr.256545.119</li>",
        "<li>Abyzov, A., Urban, A. E., Snyder, M., & Gerstein, M. (2011). CNVnator: An approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing. Genome Research, 21(6), 974–984. https://doi.org/10.1101/gr.114876.110</li>",
        "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
        "<li>Chen, X., Schulz-Trieglaff, O., Shaw, R., Barnes, B., Schlesinger, F., Källberg, M., Cox, A. J., Kruglyak, S., & Saunders, C. T. (2016). Manta: Rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32(8), 1220–1222. https://doi.org/10.1093/bioinformatics/btv710</li>",
        "<li>Yun, T., Li, H., Chang, P.-C., Lin, M. F., Carroll, A., & McLean, C. Y. (2021). Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Bioinformatics, 36(24), 5582–5589. https://doi.org/10.1093/bioinformatics/btaa1081</li>",
        params.skip_eklipse ? "" : "<li>Goudenège, D., Bris, C., Hoffmann, V., Desquiret-Dumas, V., Jardel, C., Rucheton, B., Bannwarth, S., Paquis-Flucklinger, V., Lebre, A. S., Colin, E., Amati-Bonneau, P., Bonneau, D., Reynier, P., Lenaers, G., & Procaccio, V. (2019). eKLIPse: A sensitive tool for the detection and quantification of mitochondrial DNA deletions from next-generation sequencing data. Genetics in Medicine, 21(6), 1407–1416. https://doi.org/10.1038/s41436-018-0350-8</li>",
    ]
    repeat_call_text = [
        "<li>Dolzhenko, E., Deshpande, V., Schlesinger, F., Krusche, P., Petrovski, R., Chen, S., Emig-Agius, D., Gross, A., Narzisi, G., Bowman, B., Scheffler, K., van Vugt, J. J. F. A., French, C., Sanchis-Juan, A., Ibáñez, K., Tucci, A., Lajoie, B. R., Veldink, J. H., Raymond, F. L., … Eberle, M. A. (2019). ExpansionHunter: A sequence-graph-based tool to analyze variation in short tandem repeat regions. Bioinformatics, 35(22), 4754–4756. https://doi.org/10.1093/bioinformatics/btz431</li>",
        "<li>Nilsson, D., & Magnusson, M. (2021). Moonso/stranger v0.7.1 (v0.7.1) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.4548873</li>"
    ]
    if (!params.skip_snv_annotation) {
        snv_annotation_text = [
            "<li>Rentzsch, P., Schubach, M., Shendure, J., & Kircher, M. (2021). CADD-Splice—Improving genome-wide variant effect prediction using deep learning-derived splice scores. Genome Medicine, 13(1), 31. https://doi.org/10.1186/s13073-021-00835-9</li>",
            "<li>Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Research, 47(D1), D886–D894. https://doi.org/10.1093/nar/gky1016</li>",
            "<li>Pedersen, B. S., Layer, R. M., & Quinlan, A. R. (2016). Vcfanno: Fast, flexible annotation of genetic variants. Genome Biology, 17(1), 118. https://doi.org/10.1186/s13059-016-0973-5</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>",
        ]
    }
    if (!params.skip_sv_annotation) {
        sv_annotation_text = [
            "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>",
        ]
    }
    if (!params.skip_mt_annotation) {
        mt_annotation_text = [
            "<li>Rentzsch, P., Schubach, M., Shendure, J., & Kircher, M. (2021). CADD-Splice—Improving genome-wide variant effect prediction using deep learning-derived splice scores. Genome Medicine, 13(1), 31. https://doi.org/10.1186/s13073-021-00835-9</li>",
            "<li>Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Research, 47(D1), D886–D894. https://doi.org/10.1093/nar/gky1016</li>",
            "<li>Pedersen, B. S., Layer, R. M., & Quinlan, A. R. (2016). Vcfanno: Fast, flexible annotation of genetic variants. Genome Biology, 17(1), 118. https://doi.org/10.1186/s13059-016-0973-5</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Preste, R., Clima, R., & Attimonelli, M. (2019). Human mitochondrial variant annotation with HmtNote [Preprint]. Bioinformatics. https://doi.org/10.1101/600619</li>",
            "<li>Weissensteiner, H., Pacher, D., Kloss-Brandstätter, A., Forer, L., Specht, G., Bandelt, H.-J., Kronenberg, F., Salas, A., & Schönherr, S. (2016). HaploGrep 2: Mitochondrial haplogroup classification in the era of high-throughput sequencing. Nucleic Acids Research, 44(W1), W58–W63. https://doi.org/10.1093/nar/gkw233</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>",
        ]
    }
    if (!params.skip_me_annotation) {
        me_annotation_text = [
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
        ]
    }
    qc_bam_text = [
        "<li>Broad Institute. (2023). Picard Tools. In Broad Institute, GitHub repository. http://broadinstitute.github.io/picard/</li>",
        "<li>Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: Advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292–294. https://doi.org/10.1093/bioinformatics/btv566</li>",
        "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
        "<li>Kent, W. J., Zweig, A. S., Barber, G., Hinrichs, A. S., & Karolchik, D. (2010). BigWig and BigBed: Enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204–2207. https://doi.org/10.1093/bioinformatics/btq351</li>",
        "<li>Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: Quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867–868. https://doi.org/10.1093/bioinformatics/btx699</li>",
    ]
    me_calls_text = [
        "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
        "<li>Keane, T. M., Wong, K., & Adams, D. J. (2013). RetroSeq: Transposable element discovery from next-generation sequencing data. Bioinformatics, 29(3), 389–390. https://doi.org/10.1093/bioinformatics/bts697</li>",
    ]
    preprocessing_text = [
        params.skip_fastqc ? "" : "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        params.skip_fastp  ? "" : "<li>Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta, 2(2), e107. https://doi.org/10.1002/imt2.107</li>",
    ]
    other_citation_text = [
        "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008</li>",
        "<li>McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354</li>",
        params.skip_peddy ? "" : "<li>Pedersen, B. S., & Quinlan, A. R. (2017). Who’s Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy. The American Journal of Human Genetics, 100(3), 406–413. https://doi.org/10.1016/j.ajhg.2017.01.017</li>",
        params.run_rtgvcfeval ? "<li>Cleary, J. G., Braithwaite, R., Gaastra, K., Hilbush, B. S., Inglis, S., Irvine, S. A., Jackson, A., Littin, R., Rathod, M., Ware, D., Zook, J. M., Trigg, L., & Vega, F. M. D. L. (2015). Comparing Variant Call Files for Performance Benchmarking of Next-Generation Sequencing Variant Calling Pipelines (p. 023754). bioRxiv. https://doi.org/10.1101/023754</li>" : "",
        "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352</li>",
        "<li>Chen, X., Sanchis-Juan, A., French, C. E., Connell, A. J., Delon, I., Kingsbury, Z., Chawla, A., Halpern, A. L., Taft, R. J., Bentley, D. R., Butchbach, M. E. R., Raymond, F. L., & Eberle, M. A. (2020). Spinal muscular atrophy diagnosis and carrier screening from genome sequencing data. Genetics in Medicine, 22(5), 945–953. https://doi.org/10.1038/s41436-020-0754-0</li>",
        "<li>Li, H. (2011). Tabix: Fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics, 27(5), 718–719. https://doi.org/10.1093/bioinformatics/btq671</li>",
    ]

    def concat_text = align_text +
                        variant_call_text    +
                        repeat_call_text     +
                        snv_annotation_text  +
                        sv_annotation_text   +
                        mt_annotation_text   +
                        qc_bam_text          +
                        me_calls_text        +
                        me_annotation_text   +
                        preprocessing_text   +
                        other_citation_text

    def reference_text = concat_text.unique(false) { a, b -> a <=> b } - ""
    return reference_text.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
