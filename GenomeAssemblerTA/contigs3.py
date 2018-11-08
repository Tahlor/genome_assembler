import collections, sys
import os

<<<<<<< HEAD
d = {"A": "T", "T": "A", "G": "C", "C": "G"}

def twin(text):
    out = ""
    for t in text[::-1]:
        if t in d.keys():
            out += d[t]
    return out
=======
output = open("results.txt", "w")
>>>>>>> af0a04d3263b6915edd6e2999632467d3ac05f92

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
            # fw is all possible contigs that could be connected
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break

        cand_list = [x for x in fw(c_fw[-1]) if x in d]

        # This should be 1 by construction
        assert len(cand_list)==1

        cand = cand_list[0]
        if cand == km or cand in c_fw:
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


<<<<<<< HEAD
def main_call(path, depth=0, k_range = None):
=======

def calc_n50(contigs):
    contigs.sort()
    contigs.sort(key=len, reverse=True)
    total_length = sum(len(contig) for contig in contigs)
    half_length = total_length / 2
    count = 0
    for contig in contigs:
        count += len(contig)
        if count >= half_length:
            return len(contig)


def longest_contig(contigs):
    if len(contigs) == 0:
        return 0
    return max([len(contig) for contig in contigs])


def main_call(path, depth=0):
>>>>>>> af0a04d3263b6915edd6e2999632467d3ac05f92
    """
    path: path of fasta file
    depth: kmers with frequence of less than "depth" are dropped (assumed to be errant)
    """
    reads = parse_fasta(path)
    print("\t K #CONTIGS CONTIG_LIST")
<<<<<<< HEAD
    if k_range is None:
        k_range = range(13,39,2)
    for k in k_range:
        d = build_de_bruijn(reads, k, depth)
        G, cs = all_contigs(d, k)
        print("\t", k, len(cs), cs) # cs is list of contigs
        #checker(reads, cs)
=======
    output.write("\tK\t#CONTIGS\tN50\t\tLARGEST\t\t\tCONTIG_LIST \n")
    for k in range(13,51,2):
        d = build_de_bruijn(kmers, k, depth)
        G, cs = all_contigs(d, k)
        print("\t", k, len(cs), cs) # cs is list of contigs
        output.write("\t" + str(k))
        output.write("\t" + str(len(cs)))
        output.write("\t\t\t" + str(calc_n50(cs)))
        output.write("\t\t" + str(longest_contig(cs)))
        output.write("\t\t\t\t" + str(cs) + "\n")
>>>>>>> af0a04d3263b6915edd6e2999632467d3ac05f92

def checker(reads, contigs):
    result = False
    for read in reads:
        for contig in contigs:
            result = read in contig
            if result:
                break
        if not result:
            print("CONTIGS don't include all reads")
            return

def main_loop():
    data_path = r"../data"
    print(os.getcwd())
    output.write(os.getcwd())
    files = os.listdir(data_path)
    files.sort()

    for f in files:
        if f[-6:] == ".fasta":
            print(f)
<<<<<<< HEAD
            main_call(os.path.join(data_path, f), depth=1)

def test():
    path = r"../data/6 real.error.large.fasta"
    reads = parse_fasta(path)
    k = 19
    d = build_de_bruijn(reads, k, 1)
    G, cs = all_contigs(d, k)
    print("\t", k, len(cs), cs)  # cs is list of contigs

if __name__ == "__main__":
    #main_loop()
    test()
=======
            output.write(f + "\n")
            main_call(os.path.join(data_path, f))
>>>>>>> af0a04d3263b6915edd6e2999632467d3ac05f92
