#!/bin/bash
# Sanitize files for PAML - replace long species names with short ones in both
# the tree and alignment files

# Set up paths
SP_NAME="/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Scripts/Data_Handling/Shorten_Species.sed"
PAML_IN_DIR="/Volumes/LaCie/Maize_Tandem_Evolution/Orthofinder/PAML_Comp_Tests"
SEL_OGS="/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Orthofinder/Selected_OGs_Compensated_TandemIDs.txt"
BKTFLT="/Volumes/LaCie/Maize_Tandem_Evolution/Orthofinder/25gap_FLT"
MARK_SCRIPT="/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Scripts/Data_Handling/Mark_Trees_For_CodeML.py"
TREES="/Volumes/LaCie/Maize_Tandem_Evolution/Orthofinder/Tree_25"

# First, make the PAML directory structure
~/anaconda_ete/bin/python ${MARK_SCRIPT} ${TREES} ${SEL_OGS} ${PAML_IN_DIR}

# For each selected orthogroup, modify its alignment input and tree input
# for short names
for i in $(cut -f 1 ${SEL_OGS})
do
    sed -i.bak -f ${SP_NAME} ${PAML_IN_DIR}/${i}/Null/${i}_Marked.tree
    sed -i.bak -f ${SP_NAME} ${PAML_IN_DIR}/${i}/Ha1/${i}_Marked.tree
    sed -i.bak -f ${SP_NAME} ${PAML_IN_DIR}/${i}/Ha2/${i}_Marked.tree
    sed -i.bak -f ${SP_NAME} ${PAML_IN_DIR}/${i}/Ha3/${i}_Marked.tree
    sed -f ${SP_NAME} ${BKTFLT}/${i}_BKT_25Flt.fa > ${PAML_IN_DIR}/${i}/Null/${i}_BKT.fa
    sed -f ${SP_NAME} ${BKTFLT}/${i}_BKT_25Flt.fa > ${PAML_IN_DIR}/${i}/Ha1/${i}_BKT.fa
    sed -f ${SP_NAME} ${BKTFLT}/${i}_BKT_25Flt.fa > ${PAML_IN_DIR}/${i}/Ha2/${i}_BKT.fa
    sed -f ${SP_NAME} ${BKTFLT}/${i}_BKT_25Flt.fa > ${PAML_IN_DIR}/${i}/Ha3/${i}_BKT.fa
done
