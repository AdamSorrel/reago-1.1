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
    subgraph = collapse_graph(subgraph, [], readDatabase)

    subgraph = merge_bifurcation(subgraph)
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
                            ';template_position=' + predecessor.split(';')[3].split('=')[0] + ',' + node_to_combine.split(';')[3].split('=')[1]

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