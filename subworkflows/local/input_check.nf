//
// Check input samplesheet and get read, sample, and case channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .set { sheet }

    ch_case_info = sheet.first()
                        .map { create_case_channel(it) }
    reads        = sheet.map { create_fastq_channels(it) }
    samples      = sheet.map { create_samples_channel(it) }

    emit:
    ch_case_info    // channel: [ case_id ]
    reads           // channel: [ val(meta), [ reads ] ]
    samples         // channel: [ sample_id, sex, phenotype, paternal_id, maternal_id, case_id ]

    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    //TODO: think about adding LB and PU, make sure only illumina will be used, ID can also contain a flowcell id
    meta.read_group   =     "\'@RG\\tID:"+ row.sample.split('_')[0] + "_" + row.fastq_1.split('/')[-1].split('_R1_')[0] + "_" + row.lane + "\\tPL:ILLUMINA\\tSM:"+row.sample.split('_')[0]+"\'"


    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row) {
    def sample       = [:]
    sample.id        = row.sample
    sample.gender    = row.gender
    sample.phenotype = row.phenotype
    sample.maternal  = row.maternal_id
    sample.paternal  = row.paternal_id
    sample.case_id   = row.case_id

    return sample
}

// Function to get a list of metadata (e.g. case id) for the case [ meta ]
def create_case_channel(LinkedHashMap row) {
    def case_info   = [:]
    case_info.id    = row.case_id

    return case_info
}
