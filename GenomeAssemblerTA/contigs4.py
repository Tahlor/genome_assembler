import collections, sys
import os

read_dict = collections.defaultdict(int)
duplicate_kmers = set()

def kmers(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]

def fw(km):
    """ Generate all possible nodes on forward edge
    """
    #dup = "ACGT" if km in duplicate_kmers else ""
    #print(dup)
    for x in 'ACGT':
        yield km[1:] + x

def build_de_bruijn(reads, k=31, limit=1):
    d = collections.defaultdict(int)
    kmers_seen = []
    # Count up kmer occurences
    for read in reads:
        read_dict[read]+=1
        duplicate = read_dict[read]>1
        for km in kmers(read, k):
            d[km] += 1
            if duplicate:
                duplicate_kmers.add(km)

    # Delete kmers that only appear once
    d1 = [x for x in d if d[x] <= limit]
    for x in d1:
        del d[x]
    return d

def contig_to_string(c):
    return c[0] + ''.join(x[-1] for x in c[1:])

def get_contig(d, km):
    """ Returns contig string, kmer_list that forms contig
    """
    c = forward_contig(d, km)
    return contig_to_string(c), c

def forward_contig(d, km):
    """ Given a kmer, get the contig
    """
    c_fw = [km]
    print("DUPS:",duplicate_kmers)
    while True:
        # there should be one edge to next kmer; if not, stop
            # fw is all possible contigs that could be connected
        print(c_fw)
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break

        cand_list = [x for x in fw(c_fw[-1]) if x in d]

        # This should be 1 by construction
        assert len(cand_list)==1

        cand = cand_list[0]

        if cand in duplicate_kmers:
            break

        if cand == km or cand in c_fw:
            break  # break out of cycles or mobius contigs

        c_fw.append(cand)

    return c_fw


def non_maximal_contigs(d,k):
    """ Get the max non branching path when starting at each kmer
    """
    done = set()
    r = []
    for x in d:
        if x not in done:
            contig_string, contig_kmer_list = get_contig(d, x)
            for y in contig_kmer_list:  # every sub-kmer is done
                done.add(y)
            r.append(contig_string)
    return r


def all_contigs_old(d, k):
    """ Get the contig when starting at each kmer
    """
    done = set()
    r = non_maximal_contigs(d,k)

    ## Find maximal kmers
    G = {}
    heads = {}
    for i, x in enumerate(r): # r includes suboptimal contigs
        G[i] = ([], [])
        heads[x[:k]] = (i, '+') # all of the kmer heads

    for i in G:
        x = r[i]
        for y in fw(x[-k:]): # check if a head matches
            if y in heads:
                G[i][0].append(heads[y])
    return G, r

def all_contigs(d, k):
    """ Get the contig when starting at each kmer
    """
    r = non_maximal_contigs(d,k)
    contigs=r[:]
    for possible_contig in r:
        for other_possible_contig in contigs:
            if possible_contig in other_possible_contig and possible_contig!=other_possible_contig:
                contigs.remove(possible_contig)
                break
    return contigs, r

def parse_fasta(path):
    with open(path, "r") as f:
        text = f.readlines()
    return [t[:-1] for i, t in enumerate(text) if i % 2 == 1]


def main_call(path, depth=0, k_range = None):
    """
    path: path of fasta file
    depth: kmers with frequence of less than "depth" are dropped (assumed to be errant)
    """
    kmers = parse_fasta(path)
    print("\t K #CONTIGS CONTIG_LIST")
    if k_range is None:
        k_range = range(13,39,2)
    for k in k_range:
        d = build_de_bruijn(kmers, k, depth)
        G, cs = all_contigs2(d, k)
        print("\t", k, len(cs), cs) # cs is list of contigs

if __name__ == "__main__":
    data_path = r"../data"
    print(os.getcwd())
    files = os.listdir(data_path)
    files.sort()
    for f in files:
        if f[-6:] == ".fasta":
            print(f)
            main_call(os.path.join(data_path, f))