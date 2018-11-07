cm = lambda d,ds: max([m(d,e) for e in ds if e != d])
m = lambda d,e: max([(s(d,e,o),o,e,d) for o in range(1-len(e),len(d))])
s = lambda d,e,o: sum([1 for p in range(max(0-o,0), min([len(e)-o, len(e), len(d)-o])) if e[p] == d[p+o]])
con = lambda x,o,s,c : c[0:max(0,o)] + s +  c[len(s)+o:]
a = lambda s, o : con(*cm(s, o)) if len(o) == 1 else a(con(*cm(s, o)), [ y for y in o if y != cm(s, o)[2]])
ah = lambda d : a(d[0],d[1:])


def read_fasta(filename):
    file = open(filename, "r")
    count = 0
    lines = []
    for line in file:
        if count % 2 != 0:
            lines.append(line.strip())
        count += 1
    return lines

kmers_1 = read_fasta("data/example.data.fasta")
kmers_2 = read_fasta("data/synthetic.example.noerror.small.fasta")
kmers_3 = read_fasta("data/synthetic.noerror.small.fasta")
kmers_4 = read_fasta("data/synthetic.noerror.large.fasta")
open("noerror_large_solution.txt", "w").write(ah(kmers_4))
