import nextflow.Nextflow

class CustomFunctions {


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
