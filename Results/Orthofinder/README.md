# Orthofinder Resuts
We did a lot of things with Orthofinder before we settled on the appropriate
filtering criteria for tandem duplciates. The files located in the
`Old_Filtering_Criteria` directory were generated with the old critera for
tandem duplicate flitering. These criteria were

- At least 75% pairwise similarity at the nucleotide level
- At least 50% of the sequence is aligned within clusters

If the clusters passed both of these critera, it was labeled as a "true" tandem
duplicate cluster, and if it failed at least one, it was labeled as a "false"
tandem duplicate cluster. The new criteria uses an adjusted pairwise similairty
threshold.

The metric used for the new filtering criteria is a re-scaling of pairwise
similarity based on the proportion of gaps in pairwise sequence alignment.
