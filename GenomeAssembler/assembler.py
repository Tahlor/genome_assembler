import collections
import sys
sys.setrecursionlimit(10000)


def read_fasta(filename):
    file = open(filename, "r")
    count = 0
    lines = []
    for line in file:
        if count % 2 != 0:
            lines.append(line.strip())
        count += 1
    return lines


def debrujin_graph(k, reads):
    pairs_dict = {}
    for read in reads:
        kmers = kmer_composition(k-1, read)
        kmers_list = [kmer for kmer in kmers]
        for i in range(0, len(kmers_list)):
            for j in range(0, len(kmers_list)):
                if i != j and suffix(kmers_list[i]) == prefix(kmers_list[j]):
                    if kmers_list[i] not in pairs_dict.keys():
                        pairs_dict[kmers_list[i]] = {kmers_list[j]}
                    else:
                        pairs_dict[kmers_list[i]].add(kmers_list[j])
    sorted_pairs = collections.OrderedDict(sorted(pairs_dict.items()))
    return sorted_pairs



def kmer_composition(k, dna):
    kmers = {}
    for i in range(0, len(dna)-k+1):
        kmers[dna[i:i+k]] = i
    ordered = collections.OrderedDict(sorted(kmers.items()))
    return ordered.keys()


def prefix(kmer):
    return kmer[:len(kmer)-1]


def suffix(kmer):
    return kmer[1:]


def calc_degrees(graph):
    degrees = {}
    for node, neighbors in graph.items():
        degrees[node] = (0, len(neighbors))
    for _, neighbors in graph.items():
        for node in neighbors:
            if node in degrees:
                degrees[node] = (degrees[node][0] + 1, degrees[node][1])
    return degrees


def find_start_node(degrees):
    start_node = '0'
    for node, degree in degrees.items():
        if degree[0] < degree[1]:
            start_node = node
    return start_node


def find_eulerian_path(node, graph, degrees, path):
    path += [node]

    if node not in degrees or degrees[node][1] == 0:
        return path

    while len(graph[node]) > 0:
        for temp_node in graph[node]:
            break
        graph[node].remove(temp_node)
        sub_path = find_eulerian_path(temp_node, graph, degrees, [])
        path = path[:1] + sub_path + path[1:]
    return path


def write_debruijn(graph):
    output = open("debruijn.txt", "w")
    for key in graph.keys():
        output.write(key + " -> ")
        for i in range(0, len(graph[key])):
            output.write(','.join(graph[key]))
            output.write('\n')


def contig_from_path(path):
    result = path[0]
    k = len(result)
    for i in range(1, len(path)):
        result += path[i][k-1]
    return result


reads = read_fasta("data/example.data.fasta")
graph = debrujin_graph(15, reads)
write_debruijn(graph)
degrees = calc_degrees(graph)
start_node = find_start_node(degrees)
eulerian_path = find_eulerian_path(start_node, graph, degrees, [])
# open('solution.txt', 'w').write('->'.join(eulerian_path))
result = contig_from_path(eulerian_path)
print(result)
