from collections import defaultdict

def open_file(path):
    with open(path, "rb") as f:
        text = f.readlines()[0]
    return text

def k_mer(text, n):
    """
    Returns a dict with kmers and a count
    """

    d = defaultdict(int)
    for i in range(0,len(text)-n+1):
        mer = text[i:i+n]
        d[mer] += 1
    return d

def kmers_list(dna, k, alphabetic = True):
    """
    Return a list of all kmers in alphabetic order with repeats
    """
    l = []
    for i in range(0, len(dna) - k + 1):
        mer = dna[i:i + k]
        l.append(mer)
    if alphabetic:
        l.sort()
    return l



def get_longest(d):
    #sorted_list = sorted(d.iteritems(), key=lambda (k,v): (v,k), reverse=True)
    m = max([k[1] for k in d.iteritems()])
    for item in d.iteritems():
        if item[1] == m:
            yield item[0]

if __name__=="__main__":
    path = r"D:\OneDrive\Documents\Graduate School\2018.4\CS 418\HW\HW2\rosalind_ba1b.txt"
    text = open_file(path)
    #text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    #text = "ACAACTATGCATCACTATCGGGAACTATCCT"
    d = k_mer(text, 12)
    out = " ".join(list(get_longest(d)))
    print(out)
