/*
 * Copyright (c) 2021 Maksim Masterov, SURF
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

#include "dimensions.h"

void Dimensions::findInternalIndices() {

    int start_i = 0, start_j = 0;
    int end_i = elts_loc.i, end_j = elts_loc.j;


    elts_loc_with_halo = elts_loc;

    // Check the number of elements in j-th direction
    if (decomp.getNgbPid().west != EMPTY_VAL) {
        elts_loc_with_halo.i++;
    }
    if (decomp.getNgbPid().east != EMPTY_VAL) {
        elts_loc_with_halo.i++;
    }

    // Check the number of elements in j-th direction
    if (decomp.getNgbPid().south != EMPTY_VAL) {
        elts_loc_with_halo.j++;
    }
    if (decomp.getNgbPid().north != EMPTY_VAL) {
        elts_loc_with_halo.j++;
    }

    end_i = elts_loc_with_halo.i;
    end_j = elts_loc_with_halo.j;

    // Check for the index range for the internal elements
    if (decomp.getNgbPid().west != EMPTY_VAL) {
        start_i = 1;
    }
    if (decomp.getNgbPid().south != EMPTY_VAL) {
        start_j = 1;
    }
    if (decomp.getNgbPid().east != EMPTY_VAL) {
        end_i--;
    }
    if (decomp.getNgbPid().north != EMPTY_VAL) {
        end_j--;
    }

    internal_range_i.beg = start_i;
    internal_range_i.end = end_i - 1;
    internal_range_j.beg = start_j;
    internal_range_j.end = end_j - 1;
}
