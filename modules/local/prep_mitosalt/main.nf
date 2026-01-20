process PREP_MITOSALT {
    tag "$meta.id"
    label "process_low"

    input:
    path chrsizes
    tuple val(meta2), path(genomefai)
    tuple val(meta), path(hisat2index)
    tuple val(meta4), path(mtfai)
    tuple val(meta5), path(mtfasta)
    tuple val(meta3), path(lastindex)
    val breakspan
    val breakthreshold
    val cluster_threshold
    val deletion_threshold_max
    val deletion_threshold_min
    val evalue_threshold
    val exclude
    val flank
    val hplimit
    val mito_name
    val paired_distance
    val score_threshold
    val sizelimit
    val split_distance_threshold
    val split_length
   
    output:
    path "mitosalt_config.txt", emit: msconfig

    script:
    """
    echo "hisat2 = hisat2"                                        > mitosalt_config.txt
    echo "lastal = lastal"                                        >> mitosalt_config.txt
    echo "lastsp = last-split"                                    >> mitosalt_config.txt
    echo "mfcv = maf-convert"                                     >> mitosalt_config.txt
    echo "reformat = reformat.sh"                                 >> mitosalt_config.txt
    echo "samtools = samtools"                                    >> mitosalt_config.txt
    echo "sambamba = sambamba"                                    >> mitosalt_config.txt
    echo "b2fq = bamToFastq"                                      >> mitosalt_config.txt
    echo "gcov = genomeCoverageBed"                               >> mitosalt_config.txt
    echo "intersectBed = intersectBed"                            >> mitosalt_config.txt
    echo "sortBed = sortBed"                                      >> mitosalt_config.txt
    echo "clusterBed = clusterBed"                                >> mitosalt_config.txt
    echo "randomBed = randomBed"                                  >> mitosalt_config.txt
    echo "groupBy = groupBy"                                      >> mitosalt_config.txt
    echo "bg2bw = bedGraphToBigWig"                               >> mitosalt_config.txt
    echo "hsindex = ${hisat2index}/reference"                     >> mitosalt_config.txt
    echo "faindex = ${genomefai}"                                 >> mitosalt_config.txt
    echo "lastindex = ${lastindex}/reference"                     >> mitosalt_config.txt
    echo "mtfaindex = ${mtfai}"                                   >> mitosalt_config.txt
    echo "gsize = ${chrsizes}"                                    >> mitosalt_config.txt
    echo "MT_fasta = ${mtfasta}"                                  >> mitosalt_config.txt
    echo "threads = 1"                                            >> mitosalt_config.txt
    echo "refchr = ${mito_name}"                                  >> mitosalt_config.txt
    echo "exclude = ${exclude}"                                   >> mitosalt_config.txt
    echo "orihs = 16081"                                          >> mitosalt_config.txt
    echo "orihe = 407"                                            >> mitosalt_config.txt
    echo "orils = 5730"                                           >> mitosalt_config.txt
    echo "orile = 5763"                                           >> mitosalt_config.txt
    echo "score_threshold = ${score_threshold}"                   >> mitosalt_config.txt
    echo "evalue_threshold = ${evalue_threshold}"                 >> mitosalt_config.txt
    echo "split_length = ${split_length}"                         >> mitosalt_config.txt
    echo "paired_distance = ${paired_distance}"                   >> mitosalt_config.txt
    echo "deletion_threshold_min = ${deletion_threshold_min}"     >> mitosalt_config.txt
    echo "deletion_threshold_max = ${deletion_threshold_max}"     >> mitosalt_config.txt
    echo "breakthreshold = ${breakthreshold}"                     >> mitosalt_config.txt
    echo "cluster_threshold = ${cluster_threshold}"               >> mitosalt_config.txt
    echo "breakspan = ${breakspan}"                               >> mitosalt_config.txt
    echo "sizelimit = ${sizelimit}"                               >> mitosalt_config.txt
    echo "hplimit = ${hplimit}"                                   >> mitosalt_config.txt
    echo "flank = ${flank}"                                       >> mitosalt_config.txt
    echo "split_distance_threshold = ${split_distance_threshold}" >> mitosalt_config.txt
    echo "dna = yes"                                              >> mitosalt_config.txt
    echo "enriched = no"                                          >> mitosalt_config.txt
    echo "nu_mt = yes"                                            >> mitosalt_config.txt
    echo "rmtmp = no"                                             >> mitosalt_config.txt
    echo "o_mt = yes"                                             >> mitosalt_config.txt
    echo "i_del = yes"                                            >> mitosalt_config.txt
    echo "cn_mt = yes"                                            >> mitosalt_config.txt
    """

    stub:
    """
    echo "hisat2 = hisat2" > mitosalt_config.txt
    echo "lastal = lastal" >> mitosalt_config.txt
    """

}
