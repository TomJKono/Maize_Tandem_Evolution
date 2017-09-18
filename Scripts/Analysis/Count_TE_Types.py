#!/usr/bin/env python
"""Script to take the output of the TE counts script and count the number of
times each _type_ of TE (rather than family) shows up. Takes two arguments:
    1) Tandem duplicate TE file
    2) GFF of the TE database"""

#   To take arguments
import sys

tandems = sys.argv[1]
gff = sys.argv[2]

#   Read the GFF and store the association between the TE ID and the type of
#   TE that it is.
te_id_types = {}
with open(gff, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            #   Chop out the TE ID from the metadata field
            metadata = tmp[8].split(';')
            for m in metadata:
                if m.startswith('ID='):
                    te_id = m[3:]
                    break
            #   Then, get the type of element it is. This should just be the
            #   third column.
            te_type = tmp[2]
            te_id_types[te_id] = te_type


#   Then parse through the tandem duplicates TE file, counting up the number of
#   TE types that show up
tandem_tes = {}
with open(tandems, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            #   Slice out the column that has the TE IDs that are represented
            tes = line.strip().split('\t')[8]
            if tes == 'NA':
                continue
            #   Then, for each TE, figure out its type, and keep a running
            #   total of the number of TEs we find in each type
            for t in tes.split(','):
                current_type = te_id_types[t]
                if current_type not in tandem_tes:
                    tandem_tes[current_type] = 1
                else:
                    tandem_tes[current_type] += 1

print 'TE_Type\tCount'
for t in sorted(tandem_tes):
    print t + '\t' + str(tandem_tes[t])
