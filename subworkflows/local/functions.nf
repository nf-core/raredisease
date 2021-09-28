//
// Utility functions used in subworkflow files
//
// List of functions in this file
//    * getCaseMeta()
//    * ....
//    * ....

//
// Function to get the meta info at the case level
//
def getCaseMeta(sample_meta) {
    ch_case_meta =  sample_meta
                        .first()
                        .map{
                            it ->
                                new_sample_meta = it.clone()
                                new_sample_meta.id = new_sample_meta.case_id
                                [ [ 'id':new_sample_meta.id ] ] }
    return ch_case_meta
}
