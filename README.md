# Prepare SAM alignments for RSEM expression level estimation

RSEM requires that mate pair alignments are adjacent in the SAM input, and that there are no indel, gapped or discordant alignments. This script takes a SAM alignment file as input and produces a rearranged and filtered SAM output that takes all these criteria into account. This output can then be used in `rsem-calculate-expression`. Also works with multi-mapping alignment files.

If any paired reads are present without their mate, they are discarded. The output SAM should therefore, have an even number of reads/rows.

Before using this script, ensure the SAM input file:
- is sorted by name (so all multi-mapping reads are adjacent, although not necessarily with adjacent mate pairs).
- has been filtered for mate pairs only (e.g. samtools view -f 3 sam_file). Note that this this only filters based on the FLAG id, and does not actually check if there are two reads present for each alignment. Thus, after this filter there may still be individual reads present with PAIRED FLAG id's but without their mate pair present.
- has >= 2 empty rows at end (e.g. echo "" >> sam_file)
