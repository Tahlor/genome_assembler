#import contigs4 as contigs3
import contigs4 as contigs3
import os
import main

def open_file(path):
    with open(path, "r") as f:
        text = [x.strip() for x in f.readlines()]
    return text

def run(path, depth=0, k=None):
    kmers = open_file(path)
    if k is None:
        k = len(kmers[0])
    d = contigs3.build_de_bruijn(kmers, k, depth)
    G, cs = contigs3.all_contigs(d, k)
    print("\t", k, len(cs), cs) # cs is list of contigs
    return cs

if __name__ == "__main__":
    if True:
        file_path = r"D:\PyCharm Projects\CS 418\HW\HW29 - Contigs\test2.txt"
        correct_solution_path = r"D:\PyCharm Projects\CS 418\HW\HW29 - Contigs\test2_answer.txt"
    else:
        file_path = r"D:\PyCharm Projects\CS 418\HW\HW29 - Contigs\rosalind_ba3k.txt"
        correct_solution_path = r"D:\PyCharm Projects\CS 418\HW\HW29 - Contigs\solution.txt"

    file_path = r"D:\PyCharm Projects\CS 418\HW\HW29 - Contigs\test.txt"


    contigs = run(file_path)
    if False:
        main.check_solution(contigs, correct_solution_path)
