import nextflow.Nextflow

class CustomFunctions {

    // Function to generate a pedigree file
    public static File makePed(samples, outdir) {

        def case_name  = samples[0].case_id
        def outfile  = new File(outdir +"/pipeline_info/${case_name}" + '.ped')
        outfile.text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
        def samples_list = []
        for(int i = 0; i<samples.size(); i++) {
            def sample_name =  samples[i].sample
            if (!samples_list.contains(sample_name)) {
                outfile.append('\n' + [samples[i].case_id, sample_name, samples[i].paternal, samples[i].maternal, samples[i].sex, samples[i].phenotype].join('\t'));
                samples_list.add(sample_name)
            }
        }
        return outfile
    }

    // Function to get a list of metadata (e.g. case id) for the case [ meta ]
    public static LinkedHashMap createCaseChannel(List rows) {
        def case_info    = [:]
        def probands     = []
        def upd_children = []
        def father       = ""
        def mother       = ""

        for (item in rows) {
            if (item.phenotype == 2) {
                probands.add(item.sample)
            }
            if ( (item.paternal!="0") && (item.paternal!="") && (item.maternal!="0") && (item.maternal!="") ) {
                upd_children.add(item.sample)
            }
            if ( (item.paternal!="0") && (item.paternal!="") ) {
                father = item.paternal
            }
            if ( (item.maternal!="0") && (item.maternal!="") ) {
                mother = item.maternal
            }
        }

        case_info.father       = father
        case_info.mother       = mother
        case_info.probands     = probands.unique()
        case_info.upd_children = upd_children.unique()
        case_info.id           = rows[0].case_id

        return case_info
    }

}
