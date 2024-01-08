//
// https://github.com/Midnighter/nextflow-utility-services/blob/1f84278c965255fb2bbe290b8d6c9f33f71da26f/lib/CustomChannelOperators.groovy
//
/**
 * MIT License
 *
 * Copyright (c) 2022 Moritz E. Beber
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

import groovyx.gpars.dataflow.DataflowBroadcast

/**
 * Provide a collection of custom channel operators that go beyond the nextflow default.
 *
 * @author Moritz E. Beber <https://github.com/Midnighter>
 * @author Mahesh Binzer-Panchal <https://github.com/mahesh-panchal>
 */
class CustomChannelOperators {

    /**
     * Join two channels by one or more keys from a map contained in each channel.
     *
     * The channel elements are assumed to be tuples whose size is at least two.
     * Typically, the maps to join by are in the first position of the tuples.
     * Please read https://www.nextflow.io/docs/latest/operator.html#join carefully.
     *
     * @param left The left-hand side channel in the join.
     * @param right The right-hand side channel in the join.
     * @param keys A list of strings providing the map keys to compare.
     * @param leftBy The position of the map in the left channel.
     * @param rightBy The position of the map in the right channel.
     * @param joinArgs A map of keyword arguments that is passed on to the nextflow join call.
     * @return The joined channels with the map in the original position of the left channel,
     *      followed by all elements of the right channel except for the map.
     */
    static DataflowBroadcast joinOnKeys(
            DataflowBroadcast left,
            DataflowBroadcast right,
            List<String> keys,
            Integer leftBy,
            Integer rightBy,
            Map joinArgs
    ) {
        // Extract desired keys from the left map, located at `leftBy`, and prepend them.
        DataflowBroadcast newLeft = left.map { tuple ->
            extractKeys(tuple, keys, leftBy) + tuple
        }

        // Extract desired keys from the right map, located at `rightBy`, and prepend them.
        // Also drop the map itself from the right.
        DataflowBroadcast newRight = right.map { tuple ->
            extractKeys(tuple, keys, rightBy) + removeMap(tuple, rightBy)
        }

        // Set the positions to join on explicitly.
        joinArgs.by = 0..<keys.size()

        // Apply the join channel operator to the channels and finally drop the keys used for joining tuples.
        return newLeft.join(joinArgs, newRight).map { tuple -> dropKeys(tuple, keys) }
    }

    static DataflowBroadcast joinOnKeys(
            Map joinArgs,
            DataflowBroadcast left,
            DataflowBroadcast right,
            String key,
            Integer leftBy = 0,
            Integer rightBy = 0
    ) {
        return joinOnKeys(left, right, [key], leftBy, rightBy, joinArgs)
    }

    static DataflowBroadcast joinOnKeys(
            DataflowBroadcast left,
            DataflowBroadcast right,
            String key,
            Integer leftBy = 0,
            Integer rightBy = 0
    ) {
        return joinOnKeys(left, right, [key], leftBy, rightBy, [:])
    }

    static DataflowBroadcast joinOnKeys(
            Map joinArgs,
            DataflowBroadcast left,
            DataflowBroadcast right,
            List<String> keys,
            Integer leftBy = 0,
            Integer rightBy = 0
    ) {
        return joinOnKeys(left, right, keys, leftBy, rightBy, joinArgs)
    }

    static DataflowBroadcast joinOnKeys(
            DataflowBroadcast left,
            DataflowBroadcast right,
            List<String> keys,
            Integer leftBy = 0,
            Integer rightBy = 0
    ) {
        return joinOnKeys(left, right, keys, leftBy, rightBy, [:])
    }

    /**
     * Extract values from a map at given position with given keys.
     *
     * @param tuple A tuple (`List`) channel element.
     * @param keys A list of strings denoting keys in the map.
     * @param index The position of the map in the tuple.
     * @return The list of values extracted from the map.
     */
    private static List extractKeys(List tuple, List<String> keys, Integer index) {
        return tuple[index].subMap(keys).values().toList()
    }

    /**
     * Return a new tuple without the map in it.
     *
     * @param tuple A tuple (`List`) channel element.
     * @param index The position of the map in the tuple.
     * @return A copy of the list without the map.
     */
    private static List removeMap(List tuple, Integer index) {
        return tuple[0..<index] + tuple[(index + 1)..<tuple.size()]
    }

    /**
     * Drop elements corresponding to the number of keys from the head of the list.
     *
     * @param tuple A tuple (`List`) channel element.
     * @param keys A list of strings denoting keys in the map.
     * @return The given list but without the prepended values.
     */
    private static List dropKeys(List tuple, List<String> keys) {
        return tuple.drop(keys.size())
    }

}
