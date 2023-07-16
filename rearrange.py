"""
Script to rearrange input SAM/BAM alignments for input into Rsem for expression level estimation.
The reads are arranged such that mate pairs are adjacent.
- works with multi-mapping alignments
- also note that before using this script, ensure the SAM/BAM is sorted by name (so all multi-mapping
reads are adjacent, although not necessarily with adjacent mate pairs. And ensure it's been filtered
for mate pairs only (samtools view -f 3), although note that this only filters based on the FLAG id,
and does not actually check if there are 2 reads present for each alignment. Thus, after this filter
there may still be individual reads present with PAIRED FLAG id's but without their mate pair present.

Therefore, the function of this script:
- identify both mate pairs reads and write them in adjacent rows in output.
- if any 'paired' reads are present without their mate, then exclude from output. Thus, ensure that
for each alignment (or multi-mapping alignments) there are an even number of reads/rows present.

* Note, to ensure script runs correctly be sure there are >= 2 empty rows at end of input SAM file.
"""

class SamRead():
    def __init__(self, qname: str, flag: int, rname: str, pos: int, mapq: int, cigar: str,
                 rnext: str, pnext: int, tlen: int, seq: str, qual: str, whole_read: str, tags = None):
        """
        Object for each SAM read
        :param qname: query template (read) name
        :param flag: bitwise FLAG
        :param rname: reference (fasta) sequence name
        :param pos: left-most mapping position (of mate pair range, if mate is present)
        :param mapq: mapping quality
        :param cigar: CIGAR string
        :param rnext: reference name of mate/next read ('=' if rnext is same as rname)
        :param pnext: position of mate/next read
        :param tlen: template length
        :param seq: segment (read) sequence
        :param qual: ASCII of Phred-scaled base quality +33
        :param whole_read: whole row/read as string
        :param tags: list of additional mapping data (see SAM manual)
        """
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.whole_read = whole_read

        self.tags = tags # list
        self.attr = dict()
        if self.tags:
            for tag in tags:
                tg = tag.split(':')
                if tg[2].isnumeric():
                    self.attr[tg[0]] = int(tg[2])
                else:
                    self.attr[tg[0]] = tg[2]

    def __gettag__(self, tag):
        return self.attr[tag]

    def __containstag__(self, tag):
        return tag in self.attr

    def __strand__(self):
        if (self.tlen > 0) and (self.pos < self.pnext):
            return 'R1'
        elif (self.tlen < 0) and (self.pos > self.pnext):
            return 'R2'

def count_rname(item_rname: str, list_reads):
    count = 0
    for ls in list_reads:
        if ls[2] == item_rname:
            count += 1
    return count

# Remove duplicate pairs (adjacent rows) of reads
def remove_dupl_list(list_reads):
    pairs = []
    newlist = []
    for i in range(len(list_reads)):
        if i % 2 == 0:
            pairs.append([list_reads[i], list_reads[i+1]])
    for i in reversed(range(len(pairs))):
        if pairs.count(pairs[i]) > 1:
            del pairs[i]
        elif pairs.count(pairs[i]) == 1:
            rev_pair = [pairs[i][1], pairs[i][0]]
            if pairs.count(rev_pair) >= 1:
                del pairs[i]
    for pair in pairs:
        newlist.append(pair[0])
        newlist.append(pair[1])
    return newlist

# Sorts read pairs to be adjacent within each group of reads/qnames
def sort_reads(list_reads, exceptions_out): 
    newlist = []
    uniq = set()
    for ls in list_reads:
        uniq.add(ls[2]) # determine unique rnames
    for runiq in uniq:
        runiq_list = []
        for ls in list_reads:
            if ls[2] == runiq:
                runiq_list.append(ls)
        if len(runiq_list) <= 1:
            continue
        for i in range(len(runiq_list)):
            for j in range(len(runiq_list)):
                if i != j:
                    if runiq_list[i][7] == runiq_list[j][3]: # mate pair = match pos of one read with pnext of another read
                        if runiq_list[i][8] == -runiq_list[j][8]: # match tlen values (make one negative for mate pairs)
                            newlist.append(runiq_list[i])
                            newlist.append(runiq_list[j])
    if len(newlist) % 2 == 0:
        newlist = remove_dupl_list(newlist) # remove duplicate entries
    elif len(newlist) % 2 != 0:
        exceptions_out.write('Not an even number of reads for this qname: ' + list_reads[0][0])
    return newlist

def readSamFile(samfile, outfile, exceptionsout):
    """
    Read SAM file
    :param filename: SAM file name
    :param outfile: sorted output file
    :param exceptions_out: reads that were no successfully written to outfile
    """
    generate = (r for r in open(samfile, 'r'))
    mmd = dict()
    out_file = open(outfile, 'a+')
    exceptions_out = open(exceptionsout, 'a+')
    for row in generate:
        if row.startswith('@'): # SAM header
            out_file.write(row)
        elif row == '\n': # for final empty row
            if len(mmd.keys()) == 1:
                qname_previous = list(mmd)[0] # convert mmd key to list and isolate first and only element
                if len(mmd[qname_previous]) > 1:
                    # Sort previous batch of alignments with common qnames
                    mmd[qname_previous].sort(key=lambda x: x[2])  # sort by rname
                    # Then sort such that read pairs within each group of reads/qnames are also adjacent
                    sr = sort_reads(mmd[qname_previous], exceptions_out)
                    mmd = {qname_previous: sr}
                    # Check if there is an even number of occurrences of each rname in dict
                    for ls in range(len(mmd[qname_previous])):
                        # if rname count in list is even
                        if count_rname(mmd[qname_previous][ls][2], mmd[qname_previous]) % 2 == 0:
                            out_file.write(mmd[qname_previous][ls][12] + '\n')  # write whole row
                        # if rname count in list not even
                        elif count_rname(mmd[qname_previous][ls][2], mmd[qname_previous]) % 2 != 0:
                            exceptions_out.write('Rname does not occur even number of times: '
                                                 + mmd[qname_previous][ls][0] + '\t'  # qname
                                                 + mmd[qname_previous][ls][2] + '\t'  # rname
                                                 + mmd[qname_previous][ls][5] + '\n')  # cigar
                            continue
            else:
                exceptions_out.write('Length of mmd !=1 at ' + qname_previous)
        else:
            try:
                fields = row.strip().split('\t')
                qname = fields[0]
                flag = int(fields[1])
                rname = fields[2]
                pos = int(fields[3])
                mapq = int(fields[4])
                cigar = fields[5]
                rnext = fields[6]
                pnext = int(fields[7])
                tlen = int(fields[8])
                seq = fields[9]
                qual = fields[10]
                tags = None
                if len(fields) > 11: # if tags are present
                    tags = []
                    for field in fields[11:]:
                        tags.append(field.strip('\n'))

                if rname == '*':
                    continue
                if rnext != '=':
                    continue
                if tlen == 0: # to avoid reads that may be in pairs but failed to align concordantly or discordantly
                    continue
                if any(item in cigar for item in 'SDIHN'): # avoid clipping, indels & skipped regions
                    continue

                read = [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags, row.strip('\n')]
                
                # Alternatively can initialise a SamRead object
                # read = SamRead(qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,row.strip('\n'),tags)

                if mmd == {}:
                    mmd[qname] = [read]
                elif mmd != {}:
                    if len(mmd.keys()) == 1:
                        if qname in mmd.keys():
                            mmd[qname].append(read)
                            qname_previous = qname
                        elif qname not in mmd.keys():
                            qkey = list(mmd)[0]
                            if len(mmd[qkey]) > 1:
                                # Sort previous batch of alignments with common qnames
                                mmd[qname_previous].sort(key=lambda x: x[2])  # sort by rname
                                # Then sort such that read pairs within each group of reads/qnames are also adjacent
                                sr = sort_reads(mmd[qname_previous], exceptions_out)
                                mmd = {qname_previous: sr}
                                # Check if there is an even number of occurrences of each rname in dict
                                for ls in range(len(mmd[qname_previous])):
                                    # if rname count in list is even
                                    if count_rname(mmd[qname_previous][ls][2], mmd[qname_previous]) % 2 == 0:
                                        out_file.write(mmd[qname_previous][ls][12] + '\n') # write whole row
                                    # if rname count in list not even
                                    elif count_rname(mmd[qname_previous][ls][2], mmd[qname_previous]) % 2 != 0:
                                        exceptions_out.write('Rname does not occur even number of times: '
                                                             + mmd[qname_previous][ls][0] + '\t'  # qname
                                                             + mmd[qname_previous][ls][2] + '\t'  # rname
                                                             + mmd[qname_previous][ls][5] + '\n')  # cigar
                                        continue
                            # Start new dictionary with next batch of alignments
                            mmd = {qname: [read]}
                    else:
                        exceptions_out.write('Length of mmd !=1 at ' + qname)

            except RuntimeError as e:
                raise RuntimeError(row, e)

    out_file.close()
    exceptions_out.close()
