//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .set { sheet }

    reads = sheet.map { create_fastq_channels(it) }
    sample = sheet.map { create_sample_channels(it) }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
    sample // channel: [ sample, sex, phenotype, paternal_id, maternal_id, case_id ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    // TODO: add read group to the meta map

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
def create_sample_channels(LinkedHashMap row) {
    def sample = [:]
    sample.id = row.sample
    sample.gender = row.gender
    sample.phenotype = row.phenotype
    sample.maternal = row.maternal_id
    sample.paternal = row.paternal_id
    sample.case_id = row.case_id

    return sample
}
