import collections, sys
import os

def kmers(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]

def fw(km):
    """ Generate all possible nodes on forward edge
    """
    for x in 'ACGT':
        yield km[1:] + x

def build_de_bruijn(reads, k=31, limit=1):
    d = collections.defaultdict(int)

    # Count up kmer occurences
    for read in reads:
        for km in kmers(read, k):
            d[km] += 1

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

    while True:
        # there should be one edge to next kmer; if not, stop
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break

        cand_list = [x for x in fw(c_fw[-1]) if x in d]

        # This should be 1 by construction
        assert len(cand_list)==1

        cand = cand_list[0]
        if cand == km:
            break  # break out of cycles or mobius contigs

        c_fw.append(cand)

    return c_fw


def all_contigs(d, k):
    """ Get the contig when starting at each kmer
    """
    done = set()
    r = []
    for x in d:
        if x not in done:
            contig_string, contig_kmer_list = get_contig(d, x)
            for y in contig_kmer_list: # every sub-kmer is done
                done.add(y)
            r.append(contig_string)

    ## Find maximal kmers
    G = {}
    heads = {}
    for i, x in enumerate(r):
        G[i] = ([], [])
        heads[x[:k]] = (i, '+')

    for i in G:
        x = r[i]
        for y in fw(x[-k:]):
            if y in heads:
                G[i][0].append(heads[y])
    return G, r

def parse_fasta(path):
    with open(path, "r") as f:
        text = f.readlines()
    return [t[:-1] for i, t in enumerate(text) if i % 2 == 1]


def main_call(path):
    kmers = parse_fasta(path)
    print("\t K #CONTIGS CONTIG_LIST")
    for k in range(13,39,2):
        d = build_de_bruijn(kmers, k, 0)
        G, cs = all_contigs(d, k)
        print("\t", k, len(cs), cs)

if __name__ == "__main__":
    data_path = r"../data"
    print(os.getcwd())
    files = os.listdir(data_path)
    files.sort()
    for f in files:
        if f[-6:] == ".fasta":
            print(f)
            main_call(os.path.join(data_path, f))