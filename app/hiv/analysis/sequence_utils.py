# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/11/14
content:    Support module to manipulate sequences.
'''
# Functions
def pair_generator(iterable):
    '''Generator for pairs in interleaved files, such as BAM files'''
    # Note: the last item is lost if odd
    it = iter(iterable)
    while True:
        try:
            p1 = it.next()
            p2 = it.next()
            yield (p1, p2)
        except StopIteration:
            raise


def align_muscle(*seqs, **kwargs):
    '''Global alignment of sequences via MUSCLE'''
    import subprocess as sp
    from Bio import AlignIO, SeqIO
    from Bio.Align.Applications import MuscleCommandline
    muscle_cline = MuscleCommandline(diags=True, quiet=True)
    child = sp.Popen(str(muscle_cline),
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    SeqIO.write(seqs, child.stdin, "fasta")
    child.stdin.close()
    align = AlignIO.read(child.stdout, "fasta")
    child.stderr.close()
    child.stdout.close()

    if ('sort' in kwargs) and kwargs['sort']:
        from Bio.Align import MultipleSeqAlignment as MSA
        alisort = []
        for seq in seqs:
            for row in align:
                if row.id == seq.id:
                    alisort.append(row)
                    break
        align = MSA(alisort)

    return align


def extract_mapped_pairs_subsample_open(bamfile_in, n_reads, maxreads=-1, VERBOSE=0):
    '''Extract random read pairs (pointers) from an open BAM file'''
    n_pairs_tot = get_number_reads_open(bamfile_in) // 2

    if n_pairs_tot <= n_reads:
        bamfile_in.reset()
        return pair_generator(bamfile_in)

    # Limit to the first part of the file
    if maxreads == -1:
        maxreads = n_pairs_tot
    else:
        maxreads = min(n_pairs_tot, maxreads)

    # Get the random indices of the reads to store
    ind_store = np.arange(maxreads)
    np.random.shuffle(ind_store)
    ind_store = ind_store[:n_reads]
    ind_store.sort()

    if VERBOSE >= 2:
        print 'Random indices between '+str(ind_store[0])+' and '+str(ind_store[-1]),
        print '(pairs is True)'

    output_reads = []
    n_written = 0
    for i, (read1, read2) in enumerate(pair_generator(bamfile_in)):
        if VERBOSE >= 2:
            if not ((i+1) % 10000):
                print i+1, n_written, ind_store[n_written]
    
        if i == ind_store[n_written]:
            output_reads.append((read1, read2))
            n_written += 1
    
        if n_written >= n_reads:
            break

    bamfile_in.reset()
    return output_reads

