//
// Check input samplesheet and get read, sample, and case channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow CHECK_INPUT {
    take:
        ch_samplesheet // channel: [mandatory] [ path(csv) ]

    main:
        SAMPLESHEET_CHECK ( ch_samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .set { sheet }

        case_info = sheet
                        .toList()
                        .map {create_case_channel(it)}


        reads     = sheet.map { row -> [[row.sample.split('_')[0]], row] }
                        .groupTuple()
                        .map { meta, rows ->
                            [rows, rows.size()]
                        }
                        .transpose()
                        .map { row, numLanes ->
                            create_fastq_channel(row + [num_lanes:numLanes])
                        }

        samples   = sheet.map { create_samples_channel(it) }

    emit:
        case_info       // channel: [ val(case_info) ]
        reads           // channel: [ val(meta), [ path(reads) ] ]
        samples         // channel: [ val(sample_id), val(sex), val(phenotype), val(paternal_id), val(maternal_id), val(case_id) ]
        versions  = SAMPLESHEET_CHECK.out.versions  // channel: [ path(versions.yml) ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta        = [:]
    meta.case_id    = row.case_id
    meta.sex        = row.sex
    meta.id         = row.sample
    meta.maternal   = row.maternal_id
    meta.paternal   = row.paternal_id
    meta.phenotype  = row.phenotype
    meta.single_end = row.single_end.toBoolean()
    meta.num_lanes  = row.num_lanes
    meta.read_group = "\'@RG\\tID:"+ row.fastq_1.split('/')[-1] + "\\tPL:ILLUMINA\\tSM:"+row.sample.split('_')[0]+"\'"


    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        error("ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}")
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            error("ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}")
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row) {
    def sample       = [:]
    sample.id        = row.sample
    sample.sex       = row.sex
    sample.phenotype = row.phenotype
    sample.maternal  = row.maternal_id
    sample.paternal  = row.paternal_id
    sample.case_id   = row.case_id

    return sample
}

// Function to get a list of metadata (e.g. case id) for the case [ meta ]
def create_case_channel(List rows) {
    def case_info    = [:]
    def probands     = []
    def upd_children = []
    def father       = ""
    def mother       = ""

    for (item in rows) {
        if (item.phenotype == "2") {
            probands.add(item.sample.split("_T")[0])
        }
        if ( (item.paternal_id!="0") && (item.paternal_id!="") && (item.maternal_id!="0") && (item.maternal_id!="") ) {
            upd_children.add(item.sample.split("_T")[0])
        }
        if ( (item.paternal_id!="0") && (item.paternal_id!="") ) {
            father = item.paternal_id
        }
        if ( (item.maternal_id!="0") && (item.maternal_id!="") ) {
            mother = item.maternal_id
        }
    }

    case_info.father       = father
    case_info.mother       = mother
    case_info.probands     = probands
    case_info.upd_children = upd_children
    case_info.id           = rows[0].case_id

    return case_info
}

