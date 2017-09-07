#!/usr/bin/env python
"""Generate a BED file with non-overlapping intervals that cover the B73v4
genome. Takes one argument:
    1) Inverval size
"""

import sys

# Define the chromosome sizes as a constant
CHR_SIZES = {
    '1': 307041717,
    '2': 244442276,
    '3': 235667834,
    '4': 246994605,
    '5': 223902240,
    '6': 174033170,
    '7': 182381542,
    '8': 181122637,
    '9': 159769782,
    '10': 150982314}


def build_intervals(step):
    """Return a list of intervals that cover the genome in the specified size.
    The intervals will not overlap."""
    # This is super ugly, but Python doesn't have a built-in natural sort func
    # and I don't want to write one right now.
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    intervals = []
    for c in chroms:
        i_start = 0
        end = CHR_SIZES[c]
        while i_start < end:
            i_end = i_start + step
            if i_end > end:
                i_end = end
            intervals.append((c, i_start, i_end))
            i_start += step
    return intervals


def main(size):
    """Main function."""
    try:
        i_size = int(size)
    except ValueError:
        print 'Please supply an integer argument.'
        exit(1)
    ints = build_intervals(i_size)
    # Then print them out
    for i in ints:
        bed = i[0] + '\t' + str(i[1]) + '\t' + str(i[2])
        print bed


if len(sys.argv) != 2:
    print """Generate a BED file with non-overlapping intervals that cover the B73v4
genome. Takes one argument:
    1) Inverval size"""
    exit(1)
else:
    main(sys.argv[1])
