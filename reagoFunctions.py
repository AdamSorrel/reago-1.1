import os
import os.path
import time
import networkx as nx
import operator
import shutil
import redis

# Imports global variables
#import globalVariables as g

###############################################################################
# Function definitions
###################################################################################

def write_gene(data, variables):
    # Function writing an output full genes (genes that reached the max lengh defined by user)
    gene_cnt = 1
    with open(variables.full_genes_path, "w") as fOut:
        for path, gene in data:
            fOut.write(">gene_" + str(gene_cnt) + "_len=" + str(len(gene)) + "\n")
            fOut.write(gene + "\n")
            gene_cnt += 1

def write_frag(data, variables):
    # Function writing scaffolds
    frag_cnt = 1
    with open(variables.fragments_path, "w") as fOut:
        for path, gene in data:
            fOut.write(">fragment_" + str(frag_cnt) + "_len=" + str(len(gene)) + "\n")
            fOut.write(gene + "\n")
            frag_cnt += 1


def timestamp():
    # Formating time stamp for in-line output
    return time.asctime()

def remove_bubble(G):
    while True:
        bubble_removed = False
        # Retrieving nodes (e.g ['1107.2|1611.2|2403.1'])
        all_node = G.nodes()
        # Starting a set of collapse candidates (unordered list)
        collapse_candidate = set([])

        # Looping through nodes in all_nodes
        for node in all_node:
            # In the above it would be '1107.2|1611.2|2403.1'
            # Should the node not be present, skip the loop
            if node not in G.nodes():
                continue

            successors = [s for s in G.successors(node) if G.in_degree(s) == 1 and G.out_degree(s) == 1]
            # If there are no successors (like in the example above), skip the loop
            if len(successors) < 2:
                continue
            # Dictonary of successors?
            d = {}
            # Looping through the list of successors
            for successor in successors:
                # Level 2 successors (successors of the successor)
                to_node = G.successors(successor)[0] # successor has only one successor
                # If level 2 successor isn't in a dictionary, it's added
                if to_node not in d:
                    d[to_node] = []
                # Successor is added to a node
                d[to_node].append(successor)

            # If there are more then one successor coming after a node, we have a bubble!
            for to_node in [n for n in d if len(d[n]) > 1]:
                # Merging the bubble
                new_node = merge_node(d[to_node][1:], d[to_node][0], node, G, 1,g.TIP_SIZE)
                # If there is a result of merging, toggler is changed and the loop continues from the top
                if new_node:
                    bubble_removed = True
                    # Adding new node to the collapse candidates list
                    collapse_candidate.add(new_node)

        G = collapse_graph(G, list(collapse_candidate),dat.read_db,g.read_position_db)
        # If no bubbles were removed the loop is exited
        if not bubble_removed:
            break

    return G

def remove_isolated_node(G):
    for node in G.nodes():
        # There are not in and out degrees to the node and the number of reads present is lower then 5
        # (e.g. there are 3 reads in '1107.2|1611.2|2403.1')
        # OR if the read length is shorter then 105% of the input read lenth (-l option)
        if  not G.in_degree(node) and not G.out_degree(node) and \
            (n_read_in_node(node) < 5 or len(dat.read_db[node]) < g.READ_LEN * 1.05):
            # Remove the node from the network
            G.remove_node(node)

    return G

def get_branching_aid(G_orig):
    # The function used in 'get assemblie' function.

    # The reverse is a graph with the same nodes and edges but with the directions of the edges reversed.
    G = G_orig.reverse(copy = True)
    # Dictionary
    d = {}
    starting_nodes = []

    # Looping through nodes of the reversed network
    for node_str in G.nodes():
        # Node string '3352.2|838.2|2758.1'
        d[node_str] = set([node_str])
        # If there are not in degrees node is a starting node
        if G.in_degree(node_str) == 0:
            starting_nodes.append(node_str)

    # BFS
    # Looping through the starting nodes
    for starting_node in starting_nodes:
        # Starting a list of nodes
        queue = [starting_node]
        # Looping as long as the queue is not empty
        while queue != []:
            # .pop() removes and returns the last item in the list
            front = queue.pop(0)
            # Retrieving successors of the node
            successors = G.successors(front)
            # Looping through successors
            for successor in successors:
                # Not sure what's going on here.
                d[successor] = d[successor].union(d[front])
                if successor not in queue:
                    queue.append(successor)
    return d

def confidence_increment(visited_path, next_node, future_nodes):
    # Starting an empty dictionary 'd' and a weighted_num_pair_end variable
    d, weighted_num_pair_end = {}, 0
    # Looping through a visited path
    for idx, node in enumerate(visited_path):
        # To remain explored
        for read_id in node.split("|"):
            base, end = read_id.split(".")
            d[base] = len(visited_path) - idx - 1

    for node in future_nodes:
        for read_id in node.split("|"):
            base, end = read_id.split(".")

            if base in d:
                weighted_num_pair_end += 1 * (g.CONFIDENCE_BASE ** d[base])

    return weighted_num_pair_end

def get_all_path(G, future_nodes, cur_path, paths):
    last_node = cur_path[-1]
    successors = G.successors(last_node)

    # ending node, stop recursion.
    if successors == []:
        paths.append(cur_path)
        return paths
    else:
        if len(successors) > 1:
            # Confidence increment is a custom function
            candidate = sorted([[confidence_increment(cur_path, s, future_nodes[s]), s] for s in successors])
            next_node = candidate[-1][1]
        else:
            next_node = successors[0]
        return get_all_path(G, future_nodes, cur_path + [next_node], paths)

def get_contig(path, G):
    # Retrieve contig from the read database
    contig = dat.read_db[path[0]]
    for idx in range(1, len(path)):
        prev, cur = path[idx-1], path[idx]
        seq = dat.read_db[cur]
        overlap = G[prev][cur]["overlap"]
        contig += seq[overlap:]
    return contig

def get_cm_pos(path, contig):
    min_cm_st = float("inf")
    max_cm_ed = 0

    for read_id in [r for r in"|".join(path).split("|") if g.read_position_db[r] >= 0 and (g.read_position_db[r] + g.READ_LEN <= len(contig))]:
        min_cm_st = min(min_cm_st, g.cm_pos[read_id][0])
        max_cm_ed = max(max_cm_ed, g.cm_pos[read_id][1])

    return min_cm_st, max_cm_ed

def get_assemblie(G):
    future_nodes = get_branching_aid(G)
    full_genes = []
    scaffold_candidates = []

    starting_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
    for node in starting_nodes:
        paths = get_all_path(G, future_nodes, [node], [])
        for path in paths:
            contig = get_contig(path, G,dat.read_db)
            if len(contig) >= g.FULL_LENGTH:
                if g.NEED_DEFLANK:
                    st_pos = min([g.r_pos[r][0] - g.read_position_db[r] for r in path[0].split("|")])
                    ed_pos = max([len(dat.read_db_original[r]) - g.r_pos[r][1] for r in path[-1].split("|")])
                    deflanked_contig = contig[st_pos : len(contig) - ed_pos]
                else:
                    deflanked_contig = contig
                full_genes.append([path, deflanked_contig])
            else:
                m_st, m_ed,g.read_position_db = get_cm_pos(path, contig)
                if len(contig) > 120:
                    scaffold_candidates.append([path, m_st, m_ed, contig])

    return full_genes, scaffold_candidates

def conf_connect(path_1, path_2):
    # path : list sequence ID covering this particular base (e.g. ['3352.2|838.2|2758.1|2717.2'])
    # CONFIDENCE_BASE : int (e.g. 10)

    # Initializing empty dictionary and a counter
    db, n_pe = {}, 0

    for idx, node in enumerate(path_1):
        # Enumerate returns an element of list and its position, or in other words number of the cycle
        # idx : position of element / number of cycle
        # node ; node header (e.g '3352.2|838.2|2758.1|2717.2')

        # Looping over Split header into sequence tags (e.g. ['3352.2','838.2','2758.1','2717.2'])
        for read_id in node.split("|"):
            # Splitting read id (e.g ['3352', '2'])
            base, end = read_id.split(".")
            # len(path_1) : how many sequences are overlapping this particular base
            db[base] = len(path_1) - idx - 1
            # FIXME
            # THIS ALL LOOKS A BIT STRANGE WHY WOULD LEN path_1 matter?

    # Looping through nodes in path_2
    for node in path_2:
        # Splitting each node by '|'
        for read_id in node.split("|"):
            # Splitting by . (see above)
            base, end = read_id.split(".")
            # If base exists in a database
            if base in db:
                # n_pe : Result of the function
                # CONFIDENCE_BASE is put on the power of db[base]
                n_pe += 1 * (g.CONFIDENCE_BASE ** db[base])
    return n_pe

def calculate_pairwise_segment_confidence(scaffold_candidates):
    # Function called by function scaffold

    # Scaffold candidates is a list of lists of sequences with headers, start position, stop position and nucleotide sequence
    # e.g [[['3352.2|838.2|2758.1'], 825, 957, 'CGCCTATA...']]

    # Number of candidates
    n_candidate = len(scaffold_candidates)

    # Starting an empty (filled with zeros) matrix of n x n elements, where n is n_candidate
    pairwise_confidence = [[0 for i in range(n_candidate)] for j in range(n_candidate)]
    # Looping through scaffold candidates
    # i : number of candidate (generated by the enumerate method)
    # path_1 : sequence ID
    # m_st_1 : start of sequence
    # m_ed_1 : end of sequence
    # contig_1 : nucleotide sequence
    for i, [path_1, m_st_1, m_ed_1, contig_1] in enumerate(scaffold_candidates):
        # Looping to find the second pairwise candidate
        for j, [path_2, m_st_2, m_ed_2, contig_2] in enumerate(scaffold_candidates):
            # i == j : scaffold_candidate compared with itself
            # or X : the compared region is smaller then 10 (m_ed = end position, m_st = start position)
            if i == j or min(m_ed_1, m_ed_2) - max(m_st_1, m_st_2) < 10:
                # pairwise confidence is set to ZERO
                pairwise_confidence[i][j] = 0
            else:
                # Otherwise confidence is the maximum from conf_cnnect (function above)
                pairwise_confidence[i][j] = max(conf_connect(path_1, path_2,g.CONFIDENCE_BASE), conf_connect(path_2, path_1,g.CONFIDENCE_BASE))

    return pairwise_confidence

def connect_contig(seg_1, m_st_1, m_ed_1, seg_2, m_st_2, m_ed_2):
    # Connecting contigs - used by function scaffold

    # m_st : start
    # m_ed : end

    # If entire contig is embedded within the other, return the longer of the two
    if m_st_1 >= m_st_2 and m_ed_1 <= m_ed_2 or m_st_2 >= m_st_1 and m_ed_2 <= m_ed_1:
        return seg_1 if len(seg_1) >= len(seg_2) else seg_2

    # If contig 1 starts after contig 2, swap their numbers
    if m_st_1 > m_st_2:
        seg_1, seg_2 = seg_2, seg_1

    # N_MIS : Number of mismatches is considered 80% of the shorter contig
    # FIXME hardcoded mismatch percentage
    N_MIS = int(min(len(seg_1), len(seg_2)) * 0.08)
    overlap = 0

    # Looping from end end of the shorter sequence of the two,
    # Removing the first nucleotide from the 1st segment
    # and last nucleotide from the 2nd segment in each loop.
    # Loop continues until only 10 basepairs are left (e.g [12, 11, 10])

    # Removing nucleotides from beginning and end of the two sequnces until the two match
    for i in range(min(len(seg_1), len(seg_2)), 10, -1):
        suffix = seg_1[-i:]
        prefix = seg_2[:i]

        # Setting a number of mismatches
        n_mis = 0
        # Looping over all position in the segment (i)
        for j in range(i):
            # If prefix doesn't match suffix, raise number of mismatches
            if suffix[j] != prefix[j]:
                n_mis += 1
            if n_mis > N_MIS:
                break

        if n_mis <= N_MIS:
            overlap = i

    # If there is a non zero overlap, return sequences concatenated
    if overlap > 0:
        return seg_1 + seg_2[overlap:]
    else:
        # Otherwise return them with some dots whatnot
        return seg_1 + "....." + seg_2

def scaffold(scaffold_candidates):
    # Setting a toggler and a list of putative full genes
    cont = True
    full_gene = []

    while cont:
        cont = False
        candidate_next = []
        # Calculating a pairwise segment confidence (matrix of n x n where n is number of candidates)
        pairwise_confidence = calculate_pairwise_segment_confidence(scaffold_candidates)

        used_candidate = []
        # Enumerate method returns a tuple with content (row) and a number of the loop (row_idx)
        for row_idx, row in enumerate(pairwise_confidence):
            # Skipping the loop if idx already in the used_candidate list
            if row_idx in used_candidate:
                continue

            # Retrieving the maximum key value (max_conf_val) and its id (max_conf_idx)
            max_conf_idx, max_conf_val = max(enumerate(row), key = operator.itemgetter(1))
            # Retrieving a matrix row with the max conf_val and retrieving max ids from within it.
            max_conf_idx_rev = max(enumerate(pairwise_confidence[max_conf_idx]), key = operator.itemgetter(1))[0]

            # Searching for the row that contains the highest conf_idx_rev
            # AND that one isn't equal to zero
            if row_idx == max_conf_idx_rev and max_conf_val != 0:
                # Add row with the max value in the used_candidate list
                used_candidate += [row_idx, max_conf_idx]

                # Retrieving the two scaffold candidates (from the pairwise comparison)
                candidate_1 = scaffold_candidates[row_idx]
                candidate_2 = scaffold_candidates[max_conf_idx]
                # path_1 : sequence_id (e.g '3352.2|838.2|2758.1')
                # m_st_1 : starting position (e.g. 825)
                # m_ed_1 : ending postition (e.g. 957)
                # contig_1 : sequence (e.g. "CGTTCACGG...")
                path_1, m_st_1, m_ed_1, contig_1 = candidate_1
                path_2, m_st_2, m_ed_2, contig_2 = candidate_2

                # Connecting contigs
                contig_new = connect_contig(contig_1, m_st_1, m_ed_1, contig_2, m_st_2, m_ed_2)\
                        if m_st_1 < m_st_2 else connect_contig(contig_1, m_st_1, m_ed_1, contig_2, m_st_2, m_ed_2)
                path_new = path_1 + path_2 if m_st_1 < m_st_2 else path_2 + path_1
                # Finding a sequence beginning (minimum starting position) and enging (max starting position)
                m_st_new, m_ed_new = min(m_st_1, m_st_2), max(m_ed_1, m_ed_2)

                # If the contig_lenth is smaller then FULL LENGTH (user defined), append to candidates
                if len(contig_new) < g.FULL_LENGTH:
                    candidate_next.append([path_new, m_st_new, m_ed_new, contig_new])
                else:
                # Otherwise append to full length genes
                    full_gene.append([path_new, contig_new])

                cont = True
            # If scaffold candidate doesn't contain the highest conf_idx_rev, append to candidate_next
            else:
                candidate_next.append(scaffold_candidates[row_idx])

        # Dump candidate_next in the scaffold_candidates
        if cont:
            scaffold_candidates = candidate_next
        else:
            break

    return full_gene, [[path, contig] for path, du, du, contig in scaffold_candidates]


# for testing purpose

