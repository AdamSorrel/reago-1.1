import networkx as nx
from reagoFunctions import write_fa_redis
import os
import time

def create_graph_using_rj(variables, db):
    # Wrapper function for the readjoiner
    # Accepts the dereplicated read database + graph_fn ("graph")

    graph_filename = variables.rj_dir + 'graph'

    # This generates a null graph with networkx
    # More info can be found here: https://networkx.github.io/documentation/development/reference/classes.digraph.html
    # A DiGraph stores nodes and edges with optional data, or attributes.
    G = nx.DiGraph()

    #non_dup_fn = g.rj_dir + graph_fn + ".fasta"

    ### Not saving any sequences now. Old database is already saved
    # Saving sequence database in a file using width 0, which means saving the whole sequence.
    write_fa_redis(db,variables.rj_dir+'dereplicated_database/dereplicated.fasta', 0)

    # This following functions generates files 'graph.fasta', 'graph.set.des', 'graph.set.esq', 'graph.set.rit' and 'graph.set.sds'
    # 'graph.fasta' is dereplicated file
    # 'graph.set.des' are headers of dereplicated files (most likely)

    # TODO : Implement length of insert and its size and variance
    # The 'prefilter' option removes all duplicate reads (already done by reago) and encodes them
    os.system("gt readjoiner prefilter -q -des -readset " + graph_filename + " -db " + variables.rj_dir + 'dereplicated_database/dereplicated.fasta')
    # Determines all pairs suffix-prefix matches (SPMs) -l specifies an overlap of the reads and has a strong effect on the output
    os.system("gt readjoiner overlap -memlimit 100MB -l " + str(int(variables.MIN_OVERLAP)) + " -readset " + graph_filename + "> rj.log" )
    # Converts the result into a txt file that's saved in the .edge.list output
    os.system("gt readjoiner spmtest -readset " + graph_filename + ".0 -test showlist > " + graph_filename + ".edge.list")

    start = time.time()

    print("Building networkX object.")
    # Read map dictionary
    read_map = {}
    # Count
    cnt = 0
    with open(graph_filename + ".des", encoding='windows-1250') as fSetDes:
        try:
            for line in fSetDes:
                read_map[str(cnt)] = line[:-1]
                cnt += 1
        except:
            print('Issue raised with {}'.format(line))
            quit()

    with open(graph_filename + ".edge.list", encoding="windows-1250") as fEdgeList:
        for line in fEdgeList:
            # chr(45) is a dash ('-'). Specifying the ASCII code for compatibility sake
            if chr(45) in line:
                continue
            read_1, read_2, overlap = line.split(" + ")
            read_id_1, read_id_2 = read_map[read_1], read_map[read_2]
            G.add_edge(read_id_1, read_id_2, overlap=int(overlap))

    print('It took', time.time() - start, 'seconds to build the NetworkX object.')

    return G