import time
import os
import numpy as np
import main
from collections import defaultdict
import k_mer
import pandas as pd
from main import stop
from graph import Graph
import itertools

number_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
letter_dict = {0:"A", 1:"C", 2:"G", 3:"T"}
nucleotides = np.asarray(["A", "C", "G", "T"])
df_base = pd.DataFrame(nucleotides.reshape(4,1), columns=["Base"])

def smush_kmers_together(kmer_list, overlap=None):
    if overlap is None:
        overlap = len(kmer_list[0])-1 # overlap is k-1
    #kmer_list=kmer_list[::-1]
    output = kmer_list[0]

    for kmer in kmer_list[1:]:
        output += kmer[overlap:]
    return output

def create_graph(readList, k=None):
    """
    Expects a list of DNA prefix/suffix strings, adds to graph
    """
    graph = Graph("")

    if k is None:
        k = len(readList[0])

    for read in readList:
        kmers = k_mer.kmers_list(read,k)
        for kmer in kmers:
            graph.add_node_prefix_suffix(kmer)

    outputs = []
    contigs = graph.find_contigs()
    out = []
    for contig in contigs:
        text = [t.text for t in contig]
        text = smush_kmers_together(text)
        out.append(text)
    out.sort()
    #return " ".join(out)
    return " ".join(out)

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def open_file(path):
    with open(path, "r") as f:
        text = [x.strip() for x in f.readlines()]
    return text

def parse_fasta(path):
    with open(path, "r") as f:
        text = f.readlines()
    return [t[:-1] for i,t in enumerate(text) if i % 2==1]

def main_call(path):
    ## Open the file
    kmers = parse_fasta(path)
    #for k in range(7,51):
    for k in [17]:
        out = create_graph(kmers,k=k)
        print(k, out,len(out))

if __name__ == "__main__":
    data_path = r"D:\PyCharm Projects\Genome Assembly Project\data"
    for f in os.listdir(data_path):
        #if f[-6:] == ".fasta" and f=="example.data.fasta":
        #if f[-6:] == ".fasta" and f == "synthetic.example.noerror.small.fasta":
        #if f[-6:] == ".fasta" and f == "synthetic.noerror.small.fasta":
        if f[-6:] == ".fasta" and f == "synthetic.noerror.large.fasta":
            print(f)
            main_call(os.path.join(data_path, f))