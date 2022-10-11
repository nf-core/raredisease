//
// Prepare reference genome files
//

include { UNTAR as UNTAR_VCFANNO                    } from '../../modules/nf-core/untar/main'
include { BWA_INDEX                                 } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                             } from '../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX                            } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GET_CHROM_SIZES                           } from '../../modules/local/get_chrom_sizes'
include { SENTIEON_BWAINDEX                         } from '../../modules/local/sentieon/bwamemindex'

workflow PREPARE_GENOME {
    take:
        aligner             // [mandatory] params.aligner
        bwa_index           // [mandatory] bwa_index
        bwamem2_index       // [mandatory] bwamem2_index
        fasta               // [mandatory] genome.fasta
        fai                 // [mandatory] genome.fai
        variant_catalog     // [optional] variant_catalog.json
        vcfanno_resources   // [mandatory] vcfanno resource file
        bwa_index_switch    // boolean val

    main:
        ch_fasta         = file(fasta)
        ch_versions      = Channel.empty()
        ch_bwa_index     = Channel.empty()
        ch_bwamem2_index = Channel.empty()

        // Fetch aligner index or create from scratch if required
        if (aligner == "bwamem2") {
            BWAMEM2_INDEX ( [[], ch_fasta] )
            ch_bwamem2_index  = !bwamem2_index ? BWAMEM2_INDEX.out.index : [[],file(bwamem2_index)]
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }

        if ( bwa_index_switch && !bwa_index) {
            if (aligner == "sentieon") {
                SENTIEON_BWAINDEX ( [[], ch_fasta] )
                ch_bwa_index = SENTIEON_BWAINDEX.out.index
                ch_versions = ch_versions.mix(SENTIEON_BWAINDEX.out.versions)
            } else {
                BWA_INDEX ( [[], ch_fasta] )
                ch_bwa_index = BWA_INDEX.out.index
                ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
            }
        } else if (bwa_index_switch) {
            ch_bwa_index = [[],file(bwa_index)]
        }

        if (aligner != "bwamem2" && aligner != "sentieon" ) {
            exit 1, 'Please provide a valid aligner!'
        }

        if ( fai ) {
            ch_fai = file(fai)
        } else {
            ch_fai = SAMTOOLS_FAIDX ( [[], ch_fasta] )
                        .fai
                        .collect{it[1]}
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }

        // Uncompress vcfanno resources if nothing else given
        if ( params.vcfanno_resources.endsWith('.tar.gz') ) {
            ch_vcfanno_resources = UNTAR_VCFANNO ( [[],params.vcfanno_resources] ).untar
                                    .map {
                                        id, resources ->
                                            return [resources]
                                    }
            ch_versions          = ch_versions.mix(UNTAR_VCFANNO.out.versions)
        } else {
            ch_vcfanno_resources = file(vcfanno_resources)
        }

        if ( variant_catalog && file(variant_catalog, checkIfExists:true) ) {
            ch_variant_catalog = file(variant_catalog)
        } else {
            if ( params.genome == 'GRCh38' ) {
                ch_variant_catalog = file("https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/testdata/reference/variant_catalog_grch38.json")
            } else {
                ch_variant_catalog = file("https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/testdata/reference/variant_catalog_grch37.json")
            }
        }

        ch_sequence_dict = GATK_SD ( ch_fasta ).dict
        ch_chrom_sizes   = GET_CHROM_SIZES ( ch_fai ).sizes
        ch_versions      = ch_versions.mix(GET_CHROM_SIZES.out.versions)


    emit:
        bwa_index               = ch_bwa_index              // path: bwamem2/index
        bwamem2_index           = ch_bwamem2_index          // path: bwamem2/index
        chrom_sizes             = ch_chrom_sizes            // path: chrom.sizes
        fasta                   = ch_fasta                  // path: genome.fasta
        fai                     = ch_fai                    // path: genome.fasta.fai
        sequence_dict           = ch_sequence_dict
        variant_catalog         = ch_variant_catalog        // path: variant_catalog.json
        vcfanno_resources       = ch_vcfanno_resources      // channel: [ untar'd files, ]
        versions                = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
