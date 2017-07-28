from celery import Celery
import numpy as np
from createRJGraph import create_graph_using_rj

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

    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold

    correct_sequencing_error(subgraph, readDatabase)

    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold in a reverse sense
    correct_sequencing_error_reverse(subgraph, readDatabase)

    # Dictionary of subgaph reads (already happened above)
    #subgraph_read_db = {}
    # Retrieving a node (e.g. '2403.1') from the list of a subgaph nodes (e.g. ['2403.1', '1611.2', '1107.2'])
    # for node in subgraph.nodes():
    #     Building a database of reads from each of the nodes
    #     subgraph_read_db[node] = readDatabase[node]

    # Generating a DiGraph with a readjoiner using a file 'graph'.
    subgraph = create_graph_using_rj("subgraph_temp", variables, readDatabase)
    subgraph = collapse_graph(subgraph=subgraph, candidates=[], readDatabase=readDatabase)

    # Merging bifurcation if there is less then N mistakes
    # I'm taking it out for now. Sounds a bit dodgy
    # subgraph = merge_bifurcation(subgraph=subgraph, readDatabase=readDatabase, variables=variables)

    subgraph = remove_bubble(subgraph)
    # Removes a node that's shorter then 105% of read length, or does not have any branches in/out or has fewer then 5 reads
    subgraph = remove_isolated_node(subgraph)

    # Collapsing the subgraph that has been pre-treated by the above set of functions
    subgraph = collapse_graph(subgraph, [], readDatabase)

    # Retrieving full genes (genes longer then user defined value) and partial scaffolds
    full, scaf = get_assemblie(subgraph)
    # Saving full genes
    full_genes += full
    # Saving scaffolds
    scaffold_candidates += scaf

    return full_genes,scaffold_candidates

def correct_sequencing_error(subgraph, readDatabase):
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

        # Determining the max length of sequences in the alignment
        maxLength = 0
        seqCount = 0
        for sequence, name in multipleSequenceAlighnment:
            seqCount += len(name.split('|'))
            length = len(sequence)
            if length > maxLength:
                maxLength = length

        multipleSequenceAlighnmentArray = np.empty((seqCount, maxLength), dtype='<U1')

        position = 0
        for sequence, name in multipleSequenceAlighnment:
            replicatedSeqNum = len(name.split('|'))
            sequenceArray = np.fromiter(sequence, dtype='<U1')
            if len(sequence) < maxLength:
                paddingLength = maxLength - len(sequence)
                padding = np.full((paddingLength), ' ', dtype='<U1')
                sequenceArray = np.append(sequenceArray, padding)
            sequenceArrayReplicated = np.repeat([sequenceArray], repeats=replicatedSeqNum, axis=0)
            multipleSequenceAlighnmentArray[position:position+replicatedSeqNum, ] = sequenceArrayReplicated
            position += replicatedSeqNum

        # Exchanging empty space positions for true empty values (this is to count non-zero value in each column later
        multipleSequenceAlighnmentArray[multipleSequenceAlighnmentArray == ' ',] = ''

        # Gathering alignment statistics
        nonZero = np.count_nonzero(multipleSequenceAlighnmentArray, axis=0)

        for base in ['A', 'T', 'G', 'C']:
            baseCount = sum(multipleSequenceAlighnmentArray == base)
            divergentBaseCol = multipleSequenceAlighnmentArray[:,(baseCount != nonZero) & (baseCount != 0),]
            if divergentBaseCol.size != 0:
                nDivergentBases = sum((divergentBaseCol != base) & (divergentBaseCol != ''))
                if nDivergentBases[0] == 1:
                    print('Replacing a divergent base')
                    divColsBoolean = (baseCount != nonZero) & (baseCount != 0)
                    divColPosition = np.where(divColsBoolean)
                    divRowBoolean = (divergentBaseCol != base) & (divergentBaseCol != '')
                    divRowPosition = np.where(divRowBoolean)
                    multipleSequenceAlighnmentArray[divRowPosition[0][:],divColPosition[0][0]] = base
                else:
                    print('Houston, we have a problem')

    return

def correct_sequencing_error_reverse(subgraph, readDatabase):
    # G.nodes  returns a list of nodes of the network. Should there be only one node, quitting.
    if len(subgraph.nodes()) <= 1:
        return

    # Saving a sub-dicionary with the key value of each starting node
    alignment_to_starting_node = {}
    # List of starting nodes
    starting_nodes = []

    # Set is an unordered collection of items where every element is unique
    visited = set([])
    # get starting nodes
    # Looping over nodes and determining their successors, if there are none, the node is declared a 'starting node'

    for node_str in subgraph.nodes():
        if len(subgraph.successors(node_str)) == 0:
            starting_nodes.append(node_str) # every node represents a read

    # Looping over starting nodes

    for starting_node in starting_nodes:
        # Saving a sub-dicionary with the key value of each starting node
        alignment_to_starting_node[starting_node] = {}
        alignment_to_starting_node[starting_node][starting_node], min_st_pos = 0, 0

        # Queue starting with the starting node
        queue = [starting_node] # BFS

        # Looping until the queue is empty
        # This loop effectively crawls through the network starting with a starting node and crawling through predecessors
        while queue != []:
            # Returns and removes an item from a defined position
            current = queue.pop(0)
            # Visited is defined above the loop over starting nodes
            # If starting nodes have been visited, the whole loop is skipped
            if current in visited:
                continue
            else:   # Starting nodes that have not been visited are added to the 'visited' list
                visited.add(current) # ownership of cur

            # Retreaves the list of predecessors
            predecessors = subgraph.predecessors(current)
            # Predecessors are added to the queue
            queue += predecessors

            # For each element in the list of predecessors:
            for predecessor in predecessors:
                # This retrieves the overlap value of the successor from the G network
                overlap = subgraph[predecessor][current]['overlap']
                # Getting a current read length
                readLength = predecessor.split(';')[1].split('=')[1]
                # This probably determines the position of the aligned read
                alignment_to_starting_node[starting_node][predecessor] = \
                        alignment_to_starting_node[starting_node][current] - (int(readLength) - overlap)
                # Minimum starting position (probably which read is furthers away of the end node)
                min_st_pos = min(min_st_pos, alignment_to_starting_node[starting_node][predecessor])

        # Looping over all successors from a particular starting node subracting min start position
        for predecessor in alignment_to_starting_node[starting_node]:
            alignment_to_starting_node[starting_node][predecessor] -= min_st_pos

        multipleSequenceAlighnment = []
        for predecessor in alignment_to_starting_node[starting_node]:
            multipleSequenceAlighnment.append([" " * alignment_to_starting_node[starting_node][predecessor] + readDatabase[predecessor], predecessor])

        # Determining the max length of sequences in the alignment
        maxLength = 0
        seqCount = 0
        for sequence, name in multipleSequenceAlighnment:
            seqCount += len(name.split('|'))
            length = len(sequence)
            if length > maxLength:
                maxLength = length

        multipleSequenceAlighnmentArray = np.empty((seqCount, maxLength), dtype='<U1')

        position = 0
        for sequence, name in multipleSequenceAlighnment:
            replicatedSeqNum = len(name.split('|'))
            sequenceArray = np.fromiter(sequence, dtype='<U1')
            if len(sequence) < maxLength:
                paddingLength = maxLength - len(sequence)
                padding = np.full((paddingLength), ' ', dtype='<U1')
                sequenceArray = np.append(sequenceArray, padding)
            sequenceArrayReplicated = np.repeat([sequenceArray], repeats=replicatedSeqNum, axis=0)
            multipleSequenceAlighnmentArray[position:position + replicatedSeqNum, ] = sequenceArrayReplicated
            position += replicatedSeqNum

        # Exchanging empty space positions for true empty values (this is to count non-zero value in each column later
        multipleSequenceAlighnmentArray[multipleSequenceAlighnmentArray == ' ',] = ''

        # Gathering alignment statistics
        # statArray = np.empty((5,multipleSequenceAlighnmentArray.shape[1]), dtype=np.int8)

        nonZero = np.count_nonzero(multipleSequenceAlighnmentArray, axis=0)
        # Occurences of adenine in each column
        # adenine= sum(multipleSequenceAlighnmentArray == "A")
        # Occurrences should either be equal to number of non-zero elements or be 0 themselves, if not there is a divergent base
        # nonAdenine = multipleSequenceAlighnmentArray[:,(adenine != nonZero) & (adenine != 0),]
        # sum((nonAdenine != 'A') & (nonAdenine != ''))

        for base in ['A', 'T', 'G', 'C']:
            baseCount = sum(multipleSequenceAlighnmentArray == base)
            divergentBaseCol = multipleSequenceAlighnmentArray[:, (baseCount != nonZero) & (baseCount != 0), ]
            if divergentBaseCol.size != 0:
                nDivergentBases = sum((divergentBaseCol != base) & (divergentBaseCol != ''))
                if nDivergentBases[0] == 1:
                    print('Replacing a divergent base')
                    divColsBoolean = (baseCount != nonZero) & (baseCount != 0)
                    divColPosition = np.where(divColsBoolean)
                    divRowBoolean = (divergentBaseCol != base) & (divergentBaseCol != '')
                    divRowPosition = np.where(divRowBoolean)
                    multipleSequenceAlighnmentArray[divRowPosition[0][:], divColPosition[0][0]] = base
                else:
                    print('Houston, we have a problem')
    return

def collapse_graph(subgraph, candidates, readDatabase):
    # Combining nodes in the networkx produced network
    while True:
        # Starting a list of nodes
        nodes_to_combine = []

        # If the variable 'candidates' is passed as an empty list, retrieve all nodes from a network
        if not candidates:
            all_node = subgraph.nodes()
        # Otherwise work only with nodes supplied in the 'candidates' variable
        else:
            all_node = candidates

        # Looping over nodes (either from the candidates variable or all nodes from networkx)
        for node in all_node:
            # Should both IN and OUT degrees be equal to 1 (no bifurcation?)
            if subgraph.in_degree(node) == 1 and subgraph.out_degree(subgraph.predecessors(node)[0]) == 1:
                nodes_to_combine.append(node)
                # Should candidates be supplied, collapsed nodes will be removed from the list.
                if candidates:
                    candidates.remove(node)

        # If thre are no nodes to combine, loop is exited.
        if not nodes_to_combine:
            break

        # Looping through the list of nodes to combine
        for node_to_combine in nodes_to_combine:
            # Retrieving a predecessor node
            predecessor = subgraph.predecessors(node_to_combine)[0]
            # And level 2 predecessor
            predecessors_predecessors = subgraph.predecessors(predecessor)
            # Retrieving successor
            successors = subgraph.successors(node_to_combine)

            # Header is updated using '|' separating nodes (e.g. 605.2|1822.2|979.2|637.1)
            # update graph
            combined_node = predecessor.split(';')[0] + '|' + node_to_combine.split(';')[0]
            # Retrieving the value of an overlap (number)
            overlap_to_predecessor = subgraph[predecessor][node_to_combine]['overlap']

            # Noting the new length of the combined sequences

            newLength = int(predecessor.split(';')[1].split('=')[1]) + int(node_to_combine.split(';')[1].split('=')[1]) - overlap_to_predecessor

            combined_node = str(combined_node) + ';len=' + str(newLength) + ';' + predecessor.split(';')[2] +\
                            ';template_position=' + predecessor.split(';')[3].split('=')[1]

            # Adding a combined node to the networkx graph G
            subgraph.add_node(combined_node)
            # Looping over 2.level predecessors
            for predecessors_predecessor in predecessors_predecessors:
                # overlap between predecessors and 2. level predecessors
                o = subgraph[predecessors_predecessor][predecessor]['overlap']
                # Adding an edge to the network G using combined nodes and overlap value
                subgraph.add_edge(predecessors_predecessor, combined_node, overlap = o)

            # Looping over the successors
            for successor in successors:
                # Retrieving the overlap value betwenn sucessor and the node to combine
                o = subgraph[node_to_combine][successor]['overlap']
                # Adding an edge to the network G using combined nodes and the overlap value
                subgraph.add_edge(combined_node, successor, overlap = o)

            # update sequences
            # Retrieving the offset, counting from the overlap on
            offset = len(readDatabase[predecessor]) - overlap_to_predecessor

            # Updating the read positions by the offset determined above
            readPosition = int(node_to_combine.split(';')[2].split('=')[1])
            readPosition += offset

            node_to_combine_new = node_to_combine.split(';')[0] + ';' + \
                              node_to_combine.split(';')[1] + ';' +\
                              'read_postition=' + str(readPosition) + ';' \
                              + node_to_combine.split(';')[3]

            # Updating the database entry
            readDatabase[node_to_combine_new] = readDatabase[node_to_combine]
            del readDatabase[node_to_combine]

            # Retrieving the sequence of predecessor from the read_database (dictionary of all reads)
            pred_seq = readDatabase[predecessor]
            # Retrieving the sequence of node to combine from the read database (dictionary of all reads)
            node_seq = readDatabase[node_to_combine_new]
            # Combining the two reads (predecessor + overhang of the sequence node)
            combined_seq = pred_seq + node_seq[overlap_to_predecessor:]

            # Adding the combined read to the dictionary of all reads
            readDatabase[combined_node] = combined_seq

            # clean up
            # Removing predecessor and node to commbine from the network
            subgraph.remove_node(node_to_combine)
            subgraph.remove_node(predecessor)

            # Removing the just combined nodes (predecessor and node to combine) from the dictionary of all sequences
            del readDatabase[predecessor]

            # If node to combine still is in the list of nodes to combine, it's removed
            if node_to_combine in nodes_to_combine:
                nodes_to_combine.remove(node_to_combine)
            # If predecessor is in the list of nodes to combine, it will be removed
            if predecessor in nodes_to_combine:
                nodes_to_combine.remove(predecessor)

    return subgraph

def merge_bifurcation(subgraph, readDatabase, variables):
    while True:
        merged = False
        # fork out
        # Starting a collapse_candidate as an empty set
        collapse_candidate = set([])
        # Looping over nodes in the network G
        for node in subgraph.nodes():
            # If node isn't in the network G, skip the loop
            # Not sure how this ever happens, but it does. Probably for nodes removed during the cycle.
            if node not in subgraph.nodes():
                continue

            # Retrieving a set of successors (unordered list). Can be empty.
            # In a network 1-2-3, the successor of 1 is [2], 2 is [3] and 3 is [] - empty list.
            successors = set(subgraph.successors(node))
            # If there is fewer then 2 successors, skip this loop.
            if len(successors) < 2:
                continue

            # Retrieving the potential tips of the network.
            # Tips are nodes that have no successors.
            # In a network 1-2-3, 1 and 2 have out_degree == 1 and 3 have out_degree == 0
            tip_candidates = set([s for s in successors if subgraph.out_degree(s) == 0])
            if len(tip_candidates) == 0:
                # If there are no tip candidates, skip the loop.
                # So far didn't manage to encounter this.
                continue

            # Subtracting sets will remove elements of the subtrahend from minued (yields NOT-successors)
            dst_candidates = successors - tip_candidates

            # If there are no dst candidates, retrieve dst nodes
            if len(dst_candidates) == 0:
                # Looping through the dst nodes, remove a node with maximum reads
                # n_read_in_node is a custom function splits read header by '|' and counts sequence codes
                dst_node = max([[n_read_in_node(t), t] for t in tip_candidates])[1]
                # Removing the tip candidate from the tip_candidates
                tip_candidates.remove(dst_node)
            else:
                # Looping through the dst nodes, remove a node with maximum reads
                # Custom function n_read_in_node splits read header by '|' and counts sequence codes
                dst_node = max([[n_read_in_node(d), d] for d in dst_candidates])[1]  # only one dst node with maximum dst_candidates
                dst_candidates.remove(dst_node)                                      # remove dst
                # Looping over the dst_candidates list to find if there is a tip_candidate after removing the dst_node
                dst_candidates = [d for d in dst_candidates if G.out_degree(d) == 0]
                # Merging the dst_candidates with the tip_candidates
                tip_candidates = tip_candidates.union(dst_candidates)


            # Merging the dst_nodes
            merged_node = merge_node(src_list=tip_candidates, dst=dst_node, shared=node, subgraph=subgraph, readDatabase=readDatabase, variables=variables, direction=1)

            # If the input isn't empty:
            if merged_node:
                merged = True
                # Adding a node to the collapse candidate set
                collapse_candidate.add(node)

        # Calling a custom collapse function
        subgraph = collapse_graph(subgraph=subgraph, candidates=list(collapse_candidate), readDatabase=readDatabase)

        # fork in
        # Emptying a collapse_candidate set
        collapse_candidate = set([])
        # Looping through nodes in the G network
        for node in subgraph.nodes():
            # If node have already been removed previously by this loop? Possibly?
            if node not in subgraph.nodes():
                continue
            # Retrieving node of predecessors
            predecessors = set(subgraph.predecessors(node))
            # It the node has fewer then 2 predeceessors
            if len(predecessors) < 2:
                continue
            # Retrieving tip_candidates
            tip_candidates = set([p for p in predecessors if subgraph.in_degree(p) == 0])# and G.out_degree(p) == 1])
            # If we have no tip_candidates skipping the loop
            if len(tip_candidates) == 0:
                continue

            # Subtracting sets will remove elements of the subtrahend from minued (yields NOT-successors)
            dst_candidates = predecessors - tip_candidates

            # If there are no dst candidates, retrieve dst nodes
            if len(dst_candidates) == 0:
                # Looping through the dst nodes, remove a node with maximum reads
                # n_read_in_node is a custom function splits read header by '|' and counts sequence codes
                dst_node = max([[n_read_in_node(t), t] for t in tip_candidates])[1]
                # Removing the tip candidate from the tip_candidates
                tip_candidates.remove(dst_node)
            else:
                # Looping through the dst nodes, remove a node with maximum reads
                # Custom function splits read header by '|' and counts sequence codes
                dst_node = max([[n_read_in_node(d), d] for d in dst_candidates])[1]  # only one dst node
                dst_candidates.remove(dst_node)                                      # remove dst
                # Looping over the dst_candidates list to find if there is a tip_candidate after removing the dst_node
                dst_candidates = [d for d in dst_candidates if subgraph.in_degree(d) == 0]   # and G.out_degree(d) == 1]  # only if its out-deg is 0, a node will be considered tip
                tip_candidates = tip_candidates.union(dst_candidates)

            # Merging the dst_nodes
            merged_node = merge_node(src_list=tip_candidates, dst=dst_node, shared=node, subgraph=subgraph, variables=variables,readDatabase=readDatabase,direction=-1)

            # If the input isn't empty:
            if merged_node:
                merged = True
                # Adding a node to the collapse candidate set
                collapse_candidate.add(node)

        # Calling a custom collapse function
        subgraph = collapse_graph(subgraph, list(collapse_candidate),readDatabase)

        if merged == False:
            break

    subgraph = collapse_graph(subgraph, [],readDatabase)
    return subgraph

def n_read_in_node(node):
    # Retreaving the number of reads in a given node
    read_name_list = node.split(";")[0]
    read_name_list = node.split("|")
    return len(read_name_list)

def merge_node(src_list, dst, shared, subgraph, readDatabase, variables, direction):
    # src_list (e.g. tip_candidates)
    # src_list, dst and shared look like this: '1581.2|3294.2;...'
    # Number of allowed mismatches.
    # TODO : Hardcoded mismatches
    N_MIS = 3

    # dst is a header of merged sequences (e.g. '2579.2|2295.1|1042.1')
    # retrieving a sequence belonging to the header (dst_seq) (e.g. 'GATTCA...')
    dst_seq  = readDatabase[dst]
    # retrieving overlap of the 'new' node with the node to merge with (shared)
    dst_overlap = subgraph[shared][dst]['overlap']              if direction == 1 else subgraph[dst][shared]["overlap"]
    # Sequence that's overhanging from the node to be merged with
    dst_remaining  = dst_seq[dst_overlap: ]              if direction == 1 else dst_seq[ :-dst_overlap][::-1]

    # Starting list of nodes to remove and merge
    to_remove = []
    # List of nodes (e.g. tip candidates)
    to_merge  = []


    if direction == -1:
        print('Reverse direction!')

    for src in src_list:
        # Retrieving sequences of each of the listed nodes
        try:
            src_seq  = readDatabase[src]
        except:
            print("Current src: {}".format(src))
            print("This tip has already been dealt with.")
            return None
        # Determining the overlap (integer) of the node with the node to be merged with
        src_overlap = subgraph[shared][src]['overlap']          if direction == 1 else subgraph[src][shared]['overlap']
        # Sequence that's overhanging from the node to be merged with
        src_remaining  = src_seq[src_overlap: ]          if direction == 1 else src_seq[ :-src_overlap][::-1]

        # Function counting reads in header, splitting by '|'
        # Number of reads in the src input (see above what that is)
        if n_read_in_node(src) >= 1.2 * n_read_in_node(dst):
            continue

        # Number of mismatches
        mis = 0
        # Looping over the bases of overhanging part of the new read
        seqOver = []
        for i in range(min(len(src_remaining), len(dst_remaining))):
            # If the src overhang and dst overhang don't match add a mismatch
            seqOver.append(src_remaining[i])
            if src_remaining[i] != dst_remaining[i]:
                mis += 1
                # If the number of mismatches is larger then defined above, break from the loop
                if mis > N_MIS:
                    break
        print('SeqOver : {}'.format(''.join(seqOver)))
        print('Length SeqOver : {}'.format(len(seqOver)))

        # If the number of mismatches is larger then allowed above
        if mis > N_MIS:
            # If number of reads in node is larger then the tip size
            if n_read_in_node(src) < variables.TIP_SIZE:
                to_remove.append(src)
            continue

        # Offset of the sequence to overlap
        offset = dst_overlap - src_overlap if direction == 1 else ((len(dst_seq) - dst_overlap) - (len(src_seq) - src_overlap))

        # Adding the offset to the read position (of each read separately?)
        new_read_position = str(int(src.split(';')[2].split('=')[1]) + offset)
        new_read = src.split(';')[0] + ';' + src.split(';')[1] + ';' + 'read_position=' + new_read_position + ';' + src.split(';')[3]

        # Saving the new read header with the old sequence
        readDatabase[new_read] = readDatabase.pop(src) + src_remaining

        # Appending the new sequence to the list to merge
        to_merge.append(new_read)

    # Should there be nothing to remove or merge and finish the function
    if not to_remove + to_merge:
        return None

    # Looping throught the list of nodes to remove
    for n in to_remove:
        # Removing node from the network
        subgraph.remove_node(n)

    # Merging old and new nodes
    # e.g. dst : '719.1;len=100;read_position=0;template_position=765,853'
    #      to_merge : '1758.2;len=100;read_position=4;template_position=769,857'
    if to_merge:
        dst_split = dst.split(';')
        dst_name = dst_split[0]
        dst_len = num(dst_split[1].split('=')[1])
        dst_read_position = dst_split[2].split('=')[1]
        dst_template_position =dst_split[3].split('=')[1]

        to_merge_split = to_merge.split(';')
        to_merge_name = to_merge_split[0]
        to_merge_read_position = num(to_merge_split[2].split('=')[1])
        to_merge_template_position = to_merge_split[3].split('=')[1]

        new_length =  dst_len + to_merge_read_position
        new_template_position = dst_template_position.split(',')[0] + to_merge_template_position.split(',')[1]

        new_node_header = dst_name + '|' + to_merge_name + ';' + \
                          'len=' + str(new_length) + \
                          'read_position= ' + dst_read_position + \
                          'template_position= ' + new_template_position


        # Changing the label (header) in the network
        subgraph = nx.relabel_nodes(subgraph, {dst: new_node_header}, copy = False)

        # Adding the new node in the final database (THAT HAS ALREADY BEEN DONE!)
        readDatabase[new_node_header] = readDatabase.pop(dst)

        # Looping through the list of nodes to merge
        for n in to_merge:
            # Removing a node from the network
            subgraph.remove_node(n)

        return new_node
    else:
        return dst