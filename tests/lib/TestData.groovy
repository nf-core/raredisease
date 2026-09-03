/**
 * Shared test-data builders for the local subworkflow nf-tests.
 *
 * Loaded via `libDir "tests/lib"` in nf-test.config, so the static methods are
 * available inside any test's `workflow { """ ... """ }` block, e.g.
 *
 *     input[0] = channel.of([ ${TestData.sample('ACC13778A2')}, file(...) ])
 *
 * The methods return Groovy *map-literal strings* (not Map objects) so the exact
 * text is interpolated into the generated `.nf` unchanged.
 */
class TestData {

    // The minimal 9-region GIAB trio (case 'amusingmarmoset'):
    //   ACC13778A1 = affected proband, ACC13778A2 = mother, ACC13778A3 = father
    private static final Map TRIO = [
        ACC13778A1: [sex: 1, phenotype: 2, paternal: "'ACC13778A3'", maternal: "'ACC13778A2'"],
        ACC13778A2: [sex: 2, phenotype: 1, paternal: '0',            maternal: '0'],
        ACC13778A3: [sex: 1, phenotype: 1, paternal: '0',            maternal: '0'],
    ]

    /**
     * Full alignment-style sample meta for a trio member.
     * `extra` is appended verbatim, e.g. sample('ACC13778A2', "customer_id:'cust_a2'").
     */
    static String sample(String id, String extra = null) {
        def sample_map  = TRIO[id]
        def rg = "\"'@RG\\\\tID:${id}\\\\tPL:illumina\\\\tSM:${id}'\""
        def base = "id:'${id}', sample:'${id}', single_end:false, num_lanes:1, " +
                   "read_group: ${rg}, lane:1, sex:${sample_map.sex}, phenotype:${sample_map.phenotype}, " +
                   "paternal:${sample_map.paternal}, maternal:${sample_map.maternal}, case_id:'amusingmarmoset'"
        "[${extra ? "${base}, ${extra}" : base}]"
    }
}
