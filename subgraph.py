from celery import Celery

#app = Celery('subgraph', broker = 'amqp://localhost', backend='amqp')

#@app.subgraph
def subgraph_treatment(subgraph, db, variables):

    nodeList = []
    for node in subgraph.nodes():
        nodeList.append(node)

    nodeSequences = db.hmget('read_sequence', nodeList)

    readDatabase = {}
    for node, sequence in zip(nodeList, nodeSequences):
        readDatabase[node] = sequence.decode("UTF-8")

    #filename = "subgaph" + str(num) + ".eps"
    #draw_graph(subgraph, filename)
    #num += 1

    # TODO full_genes and scaffold_candidates must be moved to REDIS
    #full_genes = []
    #scaffold_candidates = []

    # TODO : Change this - bases that reappear twice identical are NOT sequencing errors
    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold
    correct_sequencing_error(subgraph, readDatabase, variables, 5)

    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold in a reverse sense
    correct_sequencing_error_reverse(subgraph, 5)

    # Dictionary of subgaph reads
    subgraph_read_db = {}
    # Retrieving a node (e.g. '2403.1') from the list of a subgaph nodes (e.g. ['2403.1', '1611.2', '1107.2'])
    for node in subgraph.nodes():
        # Building a database of reads from each of the nodes
        subgraph_read_db[node] = g.read_db[node]

    # Generating a DiGraph with a readjoiner using a file 'graph'.
    subgraph = create_graph_using_rj("subgraph_temp")
    subgraph = collapse_graph(subgraph, [])
    subgraph = merge_bifurcation(subgraph)
    subgraph = remove_bubble(subgraph)
    # Removes a node that's shorter then 105% of read length, or does not have any branches in/out or has fewer then 5 reads
    subgraph = remove_isolated_node(subgraph)

    # Collapsing the subgraph that has been pre-treated by the above set of functions
    subgraph = collapse_graph(subgraph, [])

    # Retrieving full genes (genes longer then user defined value) and partial scaffolds
    full, scaf = get_assemblie(subgraph)
    # Saving full genes
    full_genes += full
    # Saving scaffolds
    scaffold_candidates += scaf

    return full_genes,scaffold_candidates

def correct_sequencing_error(subgraph, readDatabase, variables, ratio):
    # This sequence takes the subgraph from networkx function and a variable 'ratio'

    # G.nodes  returns a list of nodes of the network. Should there be only one node, quitting.
    if len(subgraph.nodes()) <= 1:
        return


    alignment_to_starting_node = {}
    # List of starting nodes
    starting_nodes = []

    # Set is an unordered collection of items where every element is unique
    visited = set([])
    # Looping over nodes and determining their predecessors, if there are none, the node is declared a 'starting node'
    for node_str in subgraph.nodes():
        if len(subgraph.predecessors(node_str)) == 0:
            starting_nodes.append(node_str) # every node represents a read

    # Looping over starting nodes
    for starting_node in starting_nodes:
        # Saving a sub-dicionary with the key value of each starting node is a dictionary of dictionaries
        alignment_to_starting_node[starting_node] = {}
        # Setting a value of the starting_node dictionary
        # TODO : I'm really not sure why there is an embedded with twice the same key dictionary here.
        alignment_to_starting_node[starting_node][starting_node], max_start_position = 0, 0

        # Queue starting with the starting node
        queue = [starting_node] # BFS

        # Looping until the queue is empty
        # This loop effectively crawls through the network starting with a starting node and crawling through successors
        while queue != []:
            # Returns and removes an item from a defined position (first item in the list in this case)
            current = queue.pop(0)

            # Visited is defined above the loop over starting nodes
            # If starting nodes have been visited, the whole loop is skipped
            if current in visited:
                continue
            else:   # Starting nodes that have not been visited are added to the 'visited' list
                visited.add(current)

            # Retreaves the list of successors
            successors = subgraph.successors(current)
            # Successors are added to the queque
            queue += successors

            # For each element in the list of successors:
            for successor in successors:
                # This retrieves the overlap value of the successor from the G network
                overlap = subgraph[current][successor]['overlap']
                # Getting a current read length
                readLength = successor.split(';')[1].split('=')[1]
                # This determines the position of the aligned read relative to the beginning of the starting node
                alignment_to_starting_node[starting_node][successor] \
                        = alignment_to_starting_node[starting_node][current] + int(readLength) - overlap
                # Maximum starting position (probably which read is furthers away of the starting node)
                max_start_position = max(max_start_position, alignment_to_starting_node[starting_node][successor])

        # Looping over all successors from a particular starting node
        multipleSequenceAlighnment = []
        for successor in alignment_to_starting_node[starting_node]:
            # Reads are possitionned absolutely adding empty lines before them
            multipleSequenceAlighnment.append([" " * alignment_to_starting_node[starting_node][successor] + readDatabase[successor], successor])
            # MSA would be actually interesting to print out at some point

        MSAprint = []
        for msa in multipleSequenceAlighnment:
            s = "\t".join(msa)
            MSAprint.append(s)
        with open('subraph.txt', 'w') as f:
                f.write('\n'.join(MSAprint))
        quit()
        # correcting...
        # Iterating between maximum starting position to the end of the particular read
        # TODO : Implement proper read length instead of 100 vthe
        for i in range(max_start_position + 100):
            # Gathering base statistics
            composition = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0}
            # List of reads that were already accounted for in the base pair statistics
            involved_read = []

            # Aligned read : sequence of a read
            # Read id : identificator
            for aligned_read, read_id in msa: # ____ACGATGC..ACGATC 23431.CAJ.1
                # If i isn't an empty character (introduced with the alignment) AND
                # i doesn't exceed the aligned read length
                if i < len(aligned_read) and aligned_read[i] != ' ':
                    # number of overlapping reads at a current position (+1) is added to the dictionary composition of key = current basepair position
                    composition[aligned_read[i]] += len(read_id.split("|")) + 1
                    # Adding the read_id to the list of already checked  reads
                    involved_read.append(read_id)

            # Total count of every observed base
            ttl_cnt = sum(composition[k] for k in composition)

            # Determining the dominant base
            dominant = 'X'
            # Looping over bases in the composition dict
            for base in composition:
                # Retrieving base count
                base_cnt = composition[base]
                # Should the base count be higher then a pre-set ratio, dominant base is assumed
                # 3 out of 4 : 3/(4-3+1)=3/2 => Base IS NOT dominant
                # 10 out of 11 : 10/(11-10+1)=10/2=5 => Base IS dominant
                # This implies, that if 1 in 10 bases is different, it will automatically be corrected.
                # It is questionable whether this is desirable in microbial communities where abundance of organisms can vary on a scale of many folds
                if float(base_cnt) / ((ttl_cnt - base_cnt) + 1) > ratio:
                    dominant = base
                    break

            # If there is no dominant base, the rest of the loop is skipped
            if dominant == 'X': # when no dominant base
                continue

            # Looping over all ids in involved reads list (the one that was generated by splitting read ID of dereplicated reads)
            for read_id in involved_read:
                # Retreating the sequence of the read by its ID
                orig_seq = list(dat.read_db[read_id])
                # Retreaving the currently inspected base
                cur_base = orig_seq[i - alignment_to_starting_node[starting_node][read_id]]
                # Should the ratio of a current base be lower then the ERROR_CORRECTION_THRESHOLD (user input)
                if float(composition[cur_base]) / ttl_cnt < g.ERROR_CORRECTION_THRESHOLD:
                    # The base is replaced by a dominant base
                    orig_seq[i - alignment_to_starting_node[starting_node][read_id]] = dominant
                    # Joining the read back together and returning it to the read_db
                    dat.read_db[read_id] = "".join(orig_seq)
