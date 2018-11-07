import main

class Node:
    def __init__(self, text):
        self.children = []
        self.parents = []
        self.text = text
        self.edge_visited = []
        self.edge_rexplored = []

    def add_child(self, child):
        self.children.append(child)
        self.edge_visited.append(False)
        self.edge_rexplored.append(False)
        child.parents.append(self)

    def __print__(self):
        return self.text

class Graph:
    def __init__(self,dna):
        self.all_nodes = set()
        self.all_nodes_text = set()
        self.main_string = dna

    def reset_rexplored(self):
        for node in self.all_nodes:
            node.edge_rexplored = [False] * len(node.children)

    def check_node_for_unexplored_edge(self, node, mark_as_visited=True):
        if not node.children:
            return None
        for i in range(0, len(node.children)):
            if not node.edge_visited[i]:
                if mark_as_visited:
                    node.edge_visited[i] = True
                    node.edge_rexplored[i] = True
                return node.children[i]
        return None

    def check_node_for_explored_edge(self, node, mark_as_visited=True):
        if not node.children:
            return None
        for i in range(0, len(node.children)):
            if node.edge_visited[i] and not node.edge_rexplored[i]:
                node.edge_rexplored[i] = True
                return node.children[i]
        return None

    def find_unexplored_edge(self, list_of_nodes = None, mark_as_visited=True):
        if list_of_nodes is None:
            list_of_nodes = self.all_nodes
        for node in list_of_nodes:
            unexplored_node = self.check_node_for_unexplored_edge(node, mark_as_visited=mark_as_visited)
            if not unexplored_node is None:
                return node, unexplored_node
        return None, None

    def has_unexplored_edges(self):
        for node in self.all_nodes:
            for exp_status in node.edge_visited:
                if exp_status:
                    return True
        return False

    def parse_input(self, text):
        if type(text) == type(""):
            text = text.split("\n")
        for line in text:
            a, b_list = [x.strip() for x in line.split("->")]
            b_list = b_list.split(",")
            for b in b_list:
                self.add_node(a,b)
        #self.print_graph(combine_rows=True)
        #main.stop()

    def get_node(self, text):
        if text not in self.all_nodes_text:
            new_node = Node(text)
            self.all_nodes.add(new_node)
            self.all_nodes_text.add(text)
            return new_node
        else:
            for node in self.all_nodes:
                if node.text == text:
                    return node
        print(text)
        print("ERROR")

    def add_node(self, parent_text, child_text):
        #indices = [i for i, x in enumerate(my_list) if x == "whatever"]
        parent_node = self.get_node(parent_text)
        child_node = self.get_node(child_text)
        parent_node.add_child(child_node)

    def add_node_prefix_suffix(self, text):
        if text not in self.all_nodes_text:
            new_node = Node(text)
            self.all_nodes.add(new_node)
            self.all_nodes_text.add(text)
            for node in self.all_nodes:
                # suffix = prefix
                #if node.text != new_node.text:
                if node.text[1:] == new_node.text[:-1]:
                    node.add_child(new_node)
                if node.text[:-1] == new_node.text[1:] and node.text != new_node.text:
                    new_node.add_child(node)

    def print_graph(self, combine_rows=False, alphabetize_children=True):
        out = ""
        for node in self.all_nodes:
            if node.children:
                if alphabetize_children:
                    node.children.sort(key=lambda x: x.text)
                if combine_rows:
                    out += "{} -> {}\n".format(node.text, ",".join([t.text for t in node.children]))
                else:
                    for child in node.children:
                        out+= "{} -> {}\n".format(node.text, child.text)
        return out

    def get_only_children_recursively(self, node_list, first_node = False):
        last_node = node_list[-1]
        #print(last_node.text, len(last_node.parents), len(last_node.children))
        if len(last_node.children)==1: #okay if last node has more than 1 child
            [child_node] = last_node.children
            if child_node not in node_list and len(child_node.parents)==1: # ignore cycles
                node_list = self.get_only_children_recursively(node_list+[child_node])
        return node_list

    def find_contigs(self):
        contigs = []

        # Build out non-branching path for every node
        for node in self.all_nodes:
            contig = self.get_only_children_recursively([node], first_node=True)
            contigs.append(contig)

        # Remove non-maximal ones
        final_contigs = contigs[:]
        for i,contig in enumerate(contigs):
            for j,contig2 in enumerate(contigs):
                if i==j or len(contig) > len(contig2):
                    continue
                kmers = self.smush_kmers_together([x.text for x in contig])
                kmers2 = self.smush_kmers_together([x.text for x in contig2])

                if kmers==kmers2 and len(contig)>1:
                    final_contigs.remove(contig)
                #if set(contig).issubset(contig2) and contig in final_contigs:
                # Remove if subset of another
                # or contig == contig2[:-len(contig)]
                if (contig == contig2[-len(contig):] ) and contig in final_contigs:
                    final_contigs.remove(contig)
                    break

        print(len(final_contigs))

        return final_contigs

    def smush_kmers_together(self, kmer_list, overlap=None):
        if overlap is None:
            overlap = len(kmer_list[0]) - 1  # overlap is k-1
        # kmer_list=kmer_list[::-1]
        output = kmer_list[0]

        for kmer in kmer_list[1:]:
            output += kmer[overlap:]
        return output

    def delete_nodes(self, node_list):
        for node in node_list:
            self.all_nodes.remove(node)