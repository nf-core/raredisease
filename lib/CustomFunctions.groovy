import nextflow.Nextflow

class CustomFunctions {

    // Helper function to check if a value is neither 0, "0", nor ""
    private static boolean isNonZeroNonEmpty(value) {
        return (value instanceof String && value != "" && value != "0") ||
               (value instanceof Number && value != 0)
    }

    // Function to get a list of metadata (e.g. case id) for the case [ meta ]
    public static LinkedHashMap createCaseChannel(List rows) {
        def case_info    = [:]
        def probands     = [] as Set
        def upd_children = [] as Set
        def father       = ""
        def mother       = ""

        rows.each { item ->
            if (item?.phenotype == 2) {
                probands << item.sample
            }
            if (isNonZeroNonEmpty(item?.paternal) && isNonZeroNonEmpty(item?.maternal)) {
                upd_children << item.sample
            }
            if (isNonZeroNonEmpty(item?.paternal)) {
                father = item.paternal
            }
            if (isNonZeroNonEmpty(item?.maternal)) {
                mother = item.maternal
            }
        }

        case_info.father       = father
        case_info.mother       = mother
        case_info.probands     = probands.toList()
        case_info.upd_children = upd_children.toList()
        case_info.id           = rows[0].case_id

        return case_info
    }

}
