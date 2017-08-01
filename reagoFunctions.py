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

