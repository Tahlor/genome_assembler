import collections, sys
import os
import main

output = open("results.txt", "w")

class Assembler:
    def __init__(self, reads):
        self.reads = reads

    def kmers(self, seq, k):
        for i in range(len(seq) - k + 1):
            yield seq[i:i + k]

    def fw(self, km):
        """ Generate all possible nodes on forward edge
        """
        for x in 'ACGT':
            yield km[1:] + x

    def build_de_bruijn(self, reads, k=31, limit=1):
        d = collections.defaultdict(int)

        # Count up kmer occurences
        for read in reads:
            for km in self.kmers(read, k):
                d[km] += 1

        # Delete kmers that only appear once
        d1 = [x for x in d if d[x] <= limit]
        for x in d1:
            del d[x]
        return d

    def contig_to_string(self, c):
        return c[0] + ''.join(x[-1] for x in c[1:])

    def get_contig(self, d, km):
        """ Returns contig string, kmer_list that forms contig
        """
        c = self.forward_contig(d, km)
        return self.contig_to_string(c), c


    def forward_contig(self, d, km):
        """ Given a kmer, get the contig
        """
        c_fw = [km]
        k = len(km)
        contig = km
        while True:
            # there should be one edge to next kmer; if not, stop
            possible_edges = list(self.fw(c_fw[-1]))
            #print(list(possible_edges))
            #print([x in d for x in possible_edges])
            possible_edge_ct = sum([x in d for x in possible_edges])
            #print(possible_edge_ct)
            if possible_edge_ct > 1:
                actual_edges = []
                # if len(contig) < 100:
                if True:
                    for read in self.reads:
                        tmp = contig[-33:]
                        if tmp in read:
                            idx = read.find(tmp) + len(tmp)+1
                            if idx+1 < len(read):
                                kmer=read[idx-k:idx]
                                # print(kmer, k)
                                # print(read, contig)
                                assert len(kmer)==k
                                actual_edges=[kmer]
                                break
                    if len(actual_edges)!=1:
                        break
                else:
                    break
            elif possible_edge_ct==1:
                actual_edges = [x for x in possible_edges if x in d]
                # print("DONE")
            else:
                break

            # This should be 1 by construction
            # print(actual_edges)
            assert len(actual_edges)==1

            cand = actual_edges[0]
            if cand == km or cand in c_fw:
                break  # break out of cycles or mobius contigs

            c_fw.append(cand)
            contig+=cand[-1]
        return c_fw

    def non_maximal_contigs(self, d,k):
        """ Get the max non branching path when starting at each kmer
        """
        done = set()
        r = []
        for x in d:
            if x not in done:
                contig_string, contig_kmer_list = self.get_contig(d, x)
                for y in contig_kmer_list:  # every sub-kmer is done
                    done.add(y)
                r.append(contig_string)
        return r


    def all_contigs(self, d, k):
        """ Get the contig when starting at each kmer
        """
        r = self.non_maximal_contigs(d,k)
        contigs=r[:]
        for possible_contig in r:
            for other_possible_contig in contigs:
                #print(possible_contig, other_possible_contig)
                if possible_contig in other_possible_contig and possible_contig!=other_possible_contig:
                    contigs.remove(possible_contig)
                    break
        return r, contigs

def parse_fasta(path):
    with open(path, "r") as f:
        text = f.readlines()
    return [t[:-1] for i, t in enumerate(text) if i % 2 == 1]

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
    """
    path: path of fasta file
    depth: kmers with frequence of less than "depth" are dropped (assumed to be errant)
    """
    kmers = parse_fasta(path)
    print("\t K #CONTIGS CONTIG_LIST")
    output.write("\tK\t#CONTIGS\tN50\t\tLARGEST\t\t\tCONTIG_LIST \n")
    assembler = Assembler(kmers)
    for k in range(13,100,2):
        d = assembler.build_de_bruijn(kmers, k, depth)
        G, cs = assembler.all_contigs(d, k)
        print("\t", k, str(longest_contig(cs)), len(cs), cs) # cs is list of contigs
        output.write("\t" + str(k))
        output.write("\t" + str(len(cs)))
        output.write("\t\t\t" + str(calc_n50(cs)))
        output.write("\t\t" + str(longest_contig(cs)))
        output.write("\t\t\t\t" + str(cs) + "\n")

def test():
    path = r"../data/6 real.error.large.fasta"
    reads = parse_fasta(path)
    k = 19
    d = build_de_bruijn(reads, k, 1)
    G, cs = all_contigs(d, k)
    print("\t", k, len(cs), cs)  # cs is list of contigs

if __name__ == "__main__":
    path = r"../data/6 real.error.large.fasta"
    main_call(path, depth=2)
    #test2()
