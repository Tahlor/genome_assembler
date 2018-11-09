import os
import numpy as np
import itertools
from statistics import mode

base_pairs = {0:"A", 1:"C", 2:"G", 3:"T"}

def get_fname(path):
    (parent, child) = os.path.split(path)
    fname, ext = os.path.splitext(child)
    return fname, ext, parent

def write_file(path, text, mode="w"):
    if type(text) != type(""):
        text = str(text)
    with open(path, mode) as f:
        f.write(text)

def create_arrays(text):
    text_arr = np.asarray([t for t in text])
    zero = np.zeros(len(text_arr))
    return text_arr, zero

def open_file(path):
    with open(path, "r") as f:
        text = f.readlines()
        text = [t.strip() for t in text]
    return text

def replace_array_with_map(arr, map_dict):
    if type(map_dict)==type([]):
        map_dict = dict(zip(range(0,len(map_dict)), map_dict))
    newArray = np.copy(arr).astype(object)
    for k, v in map_dict.items(): newArray[arr == k] = v
    return newArray

def parse_fasta(path):
    with open(path, "r") as f:
        text = f.readlines()
    return "".join([t[:-1] for t in text[1:]])

def generate_random_string(size, return_array = False):
    random_string = np.random.randint(0,4,size)
    if return_array:
        return base_pairs[random_string]
    else:
        return "".join(base_pairs[random_string].tolist())


number_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
letter_dict = {0:"A", 1:"C", 2:"G", 3:"T"}

def letters_to_index(kmer):
    for l in kmer:
        yield number_dict[l]

def get_max_dict(d, use_max=True):
    func = max if use_max  else min
    f_value = func(d.values())
    args = [x for x in d if d[x]==f_value]
    return args

def generate_all_kmers(k):
    return itertools.product('ACGT', repeat=k)

def normalize_by_column(np_array):
    return np_array/np_array.sum(axis=0)

def stop():
    import time
    time.sleep(1)
    Stop

def convert_to_numpy(lines_of_text, numeric=False):
    """
    Convert list of DNAs to 2D numpy array
    """

    arr=[]
    for l in lines_of_text:
        if numeric:
            number_list = [float(x) for x in l.split()]
        else:
            number_list = [x for x in l]
        arr.append(number_list)
    arr = np.asarray(arr)
    return arr

def mode_of_list(the_list, break_ties = True):
    #return mode(the_list)
    return max(set(the_list), key=the_list.count)


def check_solution(my_answer, true_answer, separator=" ", exact=False):
    mistakes = 0
    if "/" in my_answer:
        text1 = (" ".join(open_file(my_answer)).split(separator))
    else:
        text1 = my_answer
    text2 = (" ".join(open_file(true_answer)).split(separator))

    if exact:
        if text1 != text2:
            for i,t in enumerate(text1):
                if i >= len(text2):
                    mistakes += 1
                elif t != text2[i]:
                    mistakes+=1
        print(text1)
        print(text2)
    else:
        for i in text1:
            if i not in text2:
                print("{} not in correct solution".format(i))
                mistakes+=1
        for i in text2:
            if i not in text1:
                print("{} should be in my solution".format(i))
                mistakes += 1
    print("Done checking! {} mistakes".format(mistakes))


def invert_dict(d):
    return {v: k for k, v in d.items()}


rna_to_peptide = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"", "UAG":"",
    "UGU":"C", "UGC":"C", "UGA":"", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                  "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                  "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
                  "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                  "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                  "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G", }

#peptide_to_rna = invert_dict(rna_to_peptide)


def checker(reads, contigs):
    result = False
    for read in reads:
        for contig in contigs:
            result = read in contig
            if result:
                break
        if not result:
            print("CONTIGS don't include all reads")
            print(read)
