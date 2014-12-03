# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/11/14
content:    Support module to calculate local haplotypes from a bamfile.
'''
# Functions
def build_msa(haploc, VERBOSE=0, label=''):
    '''Build multiple sequence alignment from cluster of haplotypes'''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import ambiguous_dna

    seqs = [SeqRecord(Seq(seq, ambiguous_dna),
                      id=label+'count_'+str(count)+'_rank_'+str(i),
                      name=label+'count_'+str(count)+'_rank_'+str(i),
                      description='')
            for i, (seq, count) in enumerate(haploc.most_common())]

    from .sequence_utils import align_muscle
    ali = align_muscle(*seqs, sort=True)

    return ali


def cluster_haplotypes(haplo, VERBOSE=0, min_abundance=1):
    '''Cluster haplotypes (trivial for now)'''
    from collections import Counter
    haploc = Counter()
    for (seq, count) in haplo.iteritems():
        if count >= min_abundance:
            haploc[seq] = count

    return haploc


def merge_read_pair(seq1, seq2):
    '''Merge two reads of a pair, assuming the second starts later'''
    import numpy as np
    from seqanpy import align_ladder
    (score, ali1, ali2) = align_ladder(seq1, seq2, score_gapopen=-20)
    end1 = len(ali1.rstrip('-'))
    start2 = len(ali2) - len(ali2.lstrip('-'))
    overlap_ali = np.vstack([np.fromstring(a[start2: end1], 'S1')
                             for a in (ali1, ali2)])

    overlap = overlap_ali[0]
    overlap[overlap_ali[0] != overlap_ali[1]] = 'N'
    overlap = overlap.tostring()

    seq = ali1[:start2] + overlap + ali2[end1:]
    return seq


def trim_read_roi(read, start, end):
    '''Trim a single read to a region of interest'''
    seq = []
    pos_ref = read.pos
    pos_read = 0
    for (bt, bl) in read.cigar:
        # Insertions cannot end a block, and are always copied verbatim
        # We implicitely accept leading-edge insertions
        if bt == 1:
            if pos_ref >= start:
                seq.append(read.seq[pos_read: pos_read + bl])
            pos_read += bl

        # Deletions can both start and end a block, but there's nothing to copy
        elif bt == 2:
            if pos_ref + bl >= end:
                break
            pos_ref += bl

        # Matches can both start and end blocks, and need copying
        elif bt == 0:
            if pos_ref + bl > start:
                start_inblock = max(0, start - pos_ref)
                end_inblock = min(bl, end - pos_ref)
                seq.append(read.seq[pos_read + start_inblock:
                                    pos_read + end_inblock])

                # We ignore trailing-edge insertions
                if pos_ref + bl >= end:
                    break

            pos_ref += bl
            pos_read += bl

    seq = ''.join(seq)
    return seq


def get_local_haplotypes(bamfilename, start, end, VERBOSE=0, maxreads=-1):
    '''Extract reads fully covering the region, discarding insertions'''
    import sys
    import pysam

    from collections import Counter
    haplotypes = Counter()

    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        if maxreads == -1:
            from .sequence_utils import pair_generator
            reads_iter = pair_generator(bamfile)
        else:
            from .sequence_utils import extract_mapped_pairs_subsample_open
            reads_iter =  extract_mapped_pairs_subsample_open(bamfile, maxreads,
                                                              VERBOSE=VERBOSE)

        for irp, reads in enumerate(reads_iter):
            if VERBOSE >= 2:
                if not ((irp + 1) % 10000):
                    if irp + 1 != 10000:
                        sys.stdout.write("\x1b[1A\n")
                    sys.stdout.write(str(irp + 1))
                    sys.stdout.flush()

            # Sort fwd read first: this is important because with our insert
            # size we know the fwd read starts <= the rev read
            is_fwd = reads[0].is_reverse
            reads = [reads[is_fwd], reads[not is_fwd]]

            # Check for coverage of the region
            start_fwd = reads[0].pos
            end_fwd = start_fwd + sum(bl for (bt, bl) in reads[0].cigar if bt in (0, 2))
            start_rev = reads[1].pos
            end_rev = start_rev + sum(bl for (bt, bl) in reads[1].cigar if bt in (0, 2))
            overlap_len = max(0, end_fwd - start_rev)

            # Various scenarios possible
            if start_fwd > start:
                continue

            if end_rev < end:
                continue

            # No single read covers the whole region AND (the insert has a whole
            # OR a very short overlap)
            if (end_fwd < end) and (start_rev > start) and (overlap_len < 20):
                continue

            # Now the good cases
            if (start_fwd <= start) and (end_fwd >= end):
                seq = trim_read_roi(reads[0], start, end)

            elif (start_rev <= start) and (end_rev >= end):
                seq = trim_read_roi(reads[1], start, end)

            else:
                seqs = [trim_read_roi(read, start, end) for read in reads]
                seq = merge_read_pair(*seqs)

            haplotypes[seq] += 1

    return haplotypes


def get_local_haplotypes_aligned(bamfilename, start, end, label='', VERBOSE=0, maxreads=-1):
    '''Get local haplotypes, clustered and aligned'''
    haplo = get_local_haplotypes(bamfilename, start, end, VERBOSE=VERBOSE,
                                 maxreads=maxreads)
    haploc = cluster_haplotypes(haplo, VERBOSE=VERBOSE)
    if len(haploc) == 0:
        return None

    msa = build_msa(haploc, VERBOSE=VERBOSE,
                    label=label)

    return msa
