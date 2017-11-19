from celery import Celery
import numpy as np
from createRJGraph import create_graph_using_rj
import networkx as nx

######################
import inspect
######################

def isclean(subgraph):
    for node in subgraph.nodes():
        if inspect.isbuiltin(node):
            return False
    return True

def iscleanSet(inputSet):
    for node in inputSet:
        if inspect.isbuiltin(node):
            return False
    return True

# app = Celery('subgraph', broker = 'amqp://localhost', backend='amqp')

# @app.subgraph
def subgraph_treatment(subgraph, db, variables):
    # First the sequences belonging to a subgraph are retrieved from the redis database server

    nodeSequences = db.hmget('read_sequence', subgraph.nodes())

    # A read position database needs to be constructed. readDatabase contains all reads and sequences separately and
    # will be eventually updated when some sequences start getting aligned. readPositionDb contains information about
    # the position of each read within the alignment/scaffold.Read positions are saved for each read separately
    # (dereplication is not taken in account)

    # Names and sequences are coming in binary from the redis database and are recoded to unicode-8.

    readSequenceDb = {}
    readPositionDb = {}

    readSequenceDbDereplicated = {}

    # Important: Storing only dereplicated reads for readSequenceDb, un-dereplicated reads in
    # readSequenceDbDereplicated and un-dereplicated reads in the readPositionDb
    for node, sequence in zip(subgraph.nodes(), nodeSequences):
        node = node.split(';')[0]
        if '|' in node:
            for read in node.split('|'):
                readPositionDb[read] = 0
                readSequenceDbDereplicated[read] = sequence
                print("Combined read : {}".format(node))
            readSequenceDb[node] = sequence.decode("UTF-8")
        else:
            readSequenceDb[node] = sequence.decode("UTF-8")
            readPositionDb[node] = 0
            readSequenceDbDereplicated[node] = sequence


        # This is if they need to be un-dereplicated.
        # if '|' in node:
        #    for seqName in node.split('|'):
        #        readSequenceDb[seqName] = sequence.decode("UTF-8")
        #        readPositionDb[seqName] = 0
        # else:
        #    readSequenceDb[node] = sequence.decode("UTF-8")
        #    readPositionDb[node] = 0

    full_genes = []
    scaffold_candidates = []

    if not isclean(subgraph):
        print("Not clean!")

    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold

    readSequenceDb = correct_sequencing_error(subgraph, readSequenceDb)

    if not isclean(subgraph):
        print("Not clean!")


    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold in a reverse sense
    readSequenceDb = correct_sequencing_error_reverse(subgraph, readSequenceDb)

    if not isclean(subgraph):
        print("Not clean!")


    # Dictionary of subgaph reads (already happened above)
    # subgraph_read_db = {}
    # Retrieving a node (e.g. '2403.1') from the list of a subgaph nodes (e.g. ['2403.1', '1611.2', '1107.2'])
    # for node in subgraph.nodes():
    #     Building a database of reads from each of the nodes
    #     subgraph_read_db[node] = readDatabase[node]

    # Generating a DiGraph with a readjoiner using a file 'graph'.
    subgraph = create_graph_using_rj("subgraph_temp", variables, readSequenceDb)

    if not isclean(subgraph):
        print("Not clean!")


    subgraph, readSequenceDb, readPositionDb = collapse_graph(subgraph=subgraph,
                                                              candidates=[],
                                                              readSequenceDb=readSequenceDb,
                                                              readPositionDb=readPositionDb,
                                                              purge_original_nodes=False)

    if not isclean(subgraph):
        print("Not clean!")


    # Merging bifurcation if there is less then N mistakes
    # I'm taking it out for now. Sounds a bit dodgy
    subgraph, readSequenceDb, readPositionDb = merge_bifurcation(subgraph=subgraph,
                                                                 readSequenceDb=readSequenceDb,
                                                                 readPositionDb=readPositionDb,
                                                                 variables=variables)

    if not isclean(subgraph):
        print("Not clean!")

    subgraph, readSequenceDb, readPositionDb = remove_bubble(subgraph=subgraph,
                                                             readSequenceDb=readSequenceDb,
                                                             readPositionDb=readPositionDb,
                                                             variables=variables)

    if not isclean(subgraph):
        print("Not clean!")

    # Removes a node that's shorter then 105% of read length, or does not have any branches in/out or has fewer then 5 reads
    subgraph, readSequenceDb = remove_isolated_node(subgraph=subgraph,
                                                    readSequenceDb=readSequenceDb,
                                                    variables=variables)

    # Collapsing the subgraph that has been pre-treated by the above set of functions
    subgraph, readSequenceDb, readPositionDb = collapse_graph(subgraph=subgraph,
                                                              candidates=[],
                                                              readSequenceDb=readSequenceDb,
                                                              readPositionDb=readPositionDb,
                                                              purge_original_nodes=False)

    # TODO full_genes and scaffold_candidates must be moved to REDIS
    # Retrieving full genes (genes longer then user defined value) and partial scaffolds
    full, scaf, readSequenceDb = get_assemblie(subgraph=subgraph, readSequenceDb=readSequenceDb, variables=variables)
    # Saving full genes
    full_genes += full
    # Saving scaffolds
    scaffold_candidates += scaf

    return full_genes, scaffold_candidates


def correct_sequencing_error(subgraph, readSequenceDb):
    # Takes sequences from each subgraph sequence set and compares the overlaps. Should it discover a differing base in
    # the overlap, it will attempt to correct it. It should be noted that this might be undesirable for the user
    # as it would effectively erase some differences between the population members. The performance and application of
    # this procedure remains to be determined.

    # subraph.nodes returns a list of nodes of the network. Should there be only one node the method quits.
    if len(subgraph.nodes()) <= 1:
        return None, readSequenceDb

    # A starting node needs to be determined. This is a node that has no predecessors and is therefore effectively
    #  a 'beginning' of the graph.

    alignment_to_starting_node = {}
    starting_nodes = []

    # Reminder: Set is an unordered collection of items where every element is unique
    visited = set([])

    # Looping over nodes and determining their predecessors, if there are none, the node is declared a 'starting node'
    for node_str in subgraph.nodes():
        if len(subgraph.predecessors(node_str)) == 0:
            starting_nodes.append(node_str)  # every node represents a read

    # After all starting nodes have been gathered, each one of them is treated separately. The method reconstructs
    # an alignment based on supplied overlap values and then explores the overlapping regions.

    for starting_node in starting_nodes:

        # Saving a sub-dicionary with the key value of each starting node. These are essentially all the root tips of
        # the subgraph if it were a tree.
        alignment_to_starting_node[starting_node] = {}

        # Setting a nested dictionary starting with each root tip.
        alignment_to_starting_node[starting_node][starting_node], max_start_position = 0, 0

        # Queue starting with the root tip of subgraph.
        queue = [starting_node]  # BFS

        # Looping through the queue until its empty. This effectively crawls through the network starting with a
        # the root tip and crawling through the successor all the way to the top branches.
        while queue != []:
            # Returns and removes an item from a defined position (first item in the list in this case)
            current = queue.pop(0)

            # Visited: a set of already visited nodes, shared by the entire run though the subgraph. If the starting
            # node was already visited, it is skipped.
            if current in visited:
                continue
            else:  # Starting nodes that have not been visited are added to the 'visited' list
                visited.add(current)

            # The node list of successors is saved in the variable. These are the mates that are one level higher in
            #  the tree. The successors are added to the queue.
            successors = subgraph.successors(current)
            queue += successors

            # For each element in the list of successors, the overlap value with the predecessor node is retrieved
            # as well as the read length. The alignment is saved as a new entry in the previously mentioned nested
            # dictionary under its parent root node. The position is calculated as as an addition of
            # overhang (read length - overlap) to the previous node in the dictionary.
            for successor in successors:
                # This retrieves the overlap value of the successor from the G network
                overlap = subgraph[current][successor]['overlap']
                # This determines the position of the aligned read relative to the beginning of the starting node
                alignment_to_starting_node[starting_node][successor] = \
                    alignment_to_starting_node[starting_node][current] +\
                    (len(readSequenceDb[successor.split(';')[0]]) - overlap)

                # Maximum starting position (probably which read is furthers away of the starting node)
                max_start_position = max(max_start_position, alignment_to_starting_node[starting_node][successor])

        # After the read positions have been established, the method once again loops through the acquired dictionary
        # and constructs a multiple sequence alignment. An example can be found in the 'subgraph.txt' file.
        multipleSequenceAlighnment = []
        for successor in alignment_to_starting_node[starting_node]:
            # Reads are possitionned absolutely adding empty lines before them
            multipleSequenceAlighnment.append(
                [" " * alignment_to_starting_node[starting_node][successor] + readSequenceDb[successor.split(';')[0]],
                 successor])
            # MSA would be actually interesting to print out at some point

        # Determining the end of the sequence alignment (maxLength) and the count of sequences in the alignment inclusing
        # dereplicated sequences in order to construct an alignment numpy array object.
        maxLength = 0
        seqCount = 0
        for sequence, name in multipleSequenceAlighnment:
            seqCount += len(name.split('|'))
            length = len(sequence)
            if length > maxLength:
                maxLength = length

        multipleSequenceAlighnmentArray = np.empty((seqCount, maxLength), dtype='<U1')

        # A numpy array is constructed based on the data from the data gathered in the last few steps. The number of
        # columns corresponds to maxLength and the number of rows to seqCounts. The sequences are un-dereplicated upon
        # addition to the array since the errors are determined based on sequence occurrence counts.

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

        # Alignment statistics are processed for each column of the array.

        # The amount of aligned sequences is determined for each column of the array. It is expected that the beginning
        # and the end of the alignment will have 'lower coverage' much like it is described by the broken stick model.
        # Divergent bases are identified as present in columns that contain a given base (A,T,G or C), but the sum of
        # this base does not equal sum of total bases in that colum as determined by the array "totalBasesPerColumn".
        totalBasesPerColumn = np.count_nonzero(multipleSequenceAlighnmentArray, axis=0)

        for base in ['A', 'T', 'G', 'C']:
            baseCount = sum(multipleSequenceAlighnmentArray == base)
            divergentBaseColPosition = (baseCount != totalBasesPerColumn) & (baseCount != 0)

            if sum(divergentBaseColPosition) > 0:  # If there are some columns retrieved

                for columnPosition in np.where(divergentBaseColPosition == True)[
                    0]:  # Looping through all columns that are divergent
                    basesDict = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                    column = multipleSequenceAlighnmentArray[:, columnPosition]

                    for base in basesDict:  # Testing each base in the dictionary
                        count = sum(column == base)  # Base count in the given column
                        basesDict[base] = count  # Saving a base count to the base key

                    sortedBases = sorted(basesDict.items(), key=lambda kv: kv[1], reverse=True)

                    if sortedBases[1][1] / sum(basesDict.values()) < 0.1:  # #### USER DEFINED treshold!!
                        # All places that are NOT the most abundant base (position 0 in sorted list) and at the same time
                        # are not empty spaces are going to be replaced by the most abundant base.
                        multipleSequenceAlighnmentArray[
                            (multipleSequenceAlighnmentArray[:, columnPosition] != sortedBases[0][0]) & (
                                multipleSequenceAlighnmentArray[:, columnPosition] != ''), columnPosition] = \
                        sortedBases[0][0]

                        ###############################
                        # Do something
                        ###############################

    return readSequenceDb


def correct_sequencing_error_reverse(subgraph, readSequenceDb):
    # G.nodes  returns a list of nodes of the network. Should there be only one node, quitting.
    if len(subgraph.nodes()) <= 1:
        return None, readSequenceDb

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
            starting_nodes.append(node_str)  # every node represents a read

    # Looping over starting nodes

    for starting_node in starting_nodes:
        # Saving a sub-dicionary with the key value of each starting node
        alignment_to_starting_node[starting_node] = {}
        alignment_to_starting_node[starting_node][starting_node], min_start_position = 0, 0

        # Queue starting with the starting node
        queue = [starting_node]  # BFS

        # Looping until the queue is empty
        # This loop effectively crawls through the network starting with a starting node and crawling through predecessors
        while queue != []:
            # Returns and removes an item from a defined position
            current = queue.pop(0)
            # Visited is defined above the loop over starting nodes
            # If starting nodes have been visited, the whole loop is skipped
            if current in visited:
                continue
            else:  # Starting nodes that have not been visited are added to the 'visited' list
                visited.add(current)  # ownership of cur

            # Retreaves the list of predecessors
            predecessors = subgraph.predecessors(current)
            # Predecessors are added to the queue
            queue += predecessors

            # For each element in the list of predecessors:
            for predecessor in predecessors:
                # This retrieves the overlap value of the successor from the G network
                overlap = subgraph[predecessor][current]['overlap']
                # Getting a current read length


                # readLength = predecessor.split(';')[1].split('=')[1]

                # This probably determines the position of the aligned read
                alignment_to_starting_node[starting_node][predecessor] = \
                    alignment_to_starting_node[starting_node][current] - (
                    len(readSequenceDb[predecessor.split(';')[0]]) - overlap)
                # Minimum starting position (probably which read is furthers away of the end node)
                min_start_position = min(min_start_position, alignment_to_starting_node[starting_node][predecessor])

        # Looping over all successors from a particular starting node subracting min start position
        for predecessor in alignment_to_starting_node[starting_node]:
            alignment_to_starting_node[starting_node][predecessor] -= min_start_position

        multipleSequenceAlighnment = []
        for predecessor in alignment_to_starting_node[starting_node]:
            multipleSequenceAlighnment.append([" " * alignment_to_starting_node[starting_node][predecessor] +
                                               readSequenceDb[predecessor.split(';')[0]], predecessor])

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

        totalBasesPerColumn = np.count_nonzero(multipleSequenceAlighnmentArray, axis=0)
        # Occurences of adenine in each column
        # adenine= sum(multipleSequenceAlighnmentArray == "A")
        # Occurrences should either be equal to number of non-zero elements or be 0 themselves, if not there is a divergent base
        # nonAdenine = multipleSequenceAlighnmentArray[:,(adenine != nonZero) & (adenine != 0),]
        # sum((nonAdenine != 'A') & (nonAdenine != ''))

        for base in ['A', 'T', 'G', 'C']:
            baseCount = sum(multipleSequenceAlighnmentArray == base)
            divergentBaseColPosition = (baseCount != totalBasesPerColumn) & (baseCount != 0)

            if sum(divergentBaseColPosition) > 0:  # If there are some columns retrieved

                for columnPosition in np.where(divergentBaseColPosition == True)[0]:
                    # Looping through all columns that are divergent
                    basesDict = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                    column = multipleSequenceAlighnmentArray[:, columnPosition]

                    for base in basesDict:  # Testing each base in the dictionary
                        count = sum(column == base)  # Base count in the given column
                        basesDict[base] = count  # Saving a base count to the base key

                    sortedBases = sorted(basesDict.items(), key=lambda kv: kv[1], reverse=True)

                    if sortedBases[1][1] / sum(basesDict.values()) < 0.1:  # #### USER DEFINED treshold!!
                        # All places that are NOT the most abundant base (position 0 in sorted list) and at the same time
                        # are not empty spaces are going to be replaced by the most abundant base.
                        multipleSequenceAlighnmentArray[
                            (multipleSequenceAlighnmentArray[:, columnPosition] != sortedBases[0][0]) & (
                                multipleSequenceAlighnmentArray[:, columnPosition] != ''), columnPosition] = \
                        sortedBases[0][0]


                        # if divergentBaseCol.size != 0:
                        #     nDivergentBases = sum((divergentBaseCol != base) & (divergentBaseCol != ''))
                        #     if nDivergentBases[0] == 1:
                        #         print('Replacing a divergent base')
                        #         divColsBoolean = (baseCount != nonZero) & (baseCount != 0)
                        #         divColPosition = np.where(divColsBoolean)
                        #         divRowBoolean = (divergentBaseCol != base) & (divergentBaseCol != '')
                        #         divRowPosition = np.where(divRowBoolean)
                        #         multipleSequenceAlighnmentArray[divRowPosition[0][:], divColPosition[0][0]] = base
                        #     else:
                        #         print('Houston, we have a problem')
    return readSequenceDb


def collapse_graph(subgraph, candidates, readSequenceDb, readPositionDb, purge_original_nodes):
    # Combining overlapping nodes in the networkx produced network. Nodes will not be combined with divergent bases.


    while True:
        # Starting a list of nodes to be combined.
        nodes_to_combine = []

        # If the variable 'candidates' is passed as an empty list, retrieve all nodes from a network.
        # Otherwise work only with nodes supplied in the 'candidates' variable.
        if not candidates:
            all_nodes = subgraph.nodes()
        else:
            all_nodes = candidates

        # Looping over nodes (either from the candidates variable or all nodes from networkx).
        # The IN and OUT degrees must be equal to 1, which means there is no bifurcation and the node is
        # not a starting or finishing node.
        # Should candidate nodes be supplied, the merged node is removed from the candidate list.
        for node in all_nodes:
            if subgraph.in_degree(node) == 1 and subgraph.out_degree(subgraph.predecessors(node)[0]) == 1:
                # Checking if there are no divergent bases between the node and its predecessor
                if check_divergent_bases(predecessor=subgraph.predecessors(node)[0],
                                         node=node,
                                         graph=subgraph,
                                         readSequenceDb=readSequenceDb):
                    nodes_to_combine.append(node)
                if candidates:
                    candidates.remove(node)

        # If thre are no nodes to combine, loop is interrupted.
        if not nodes_to_combine:
            break

        # Looping through the list of nodes to combine and retrieving the node surrounding; (I.) node's predecessor,
        # (II.) the predecessors' predecessor (2nd degree predecessor) and (III.) the node's successor.
        for node_to_combine in nodes_to_combine:
            print(nodes_to_combine)
            predecessor = subgraph.predecessors(node_to_combine)[0]
            predecessors_predecessors = subgraph.predecessors(predecessor)
            successors = subgraph.successors(node_to_combine)

            # If a divergent base is detected, the loop iteration is skipped
            if not check_divergent_bases(predecessor=predecessor,
                                         node=node_to_combine,
                                         graph=subgraph,
                                         readSequenceDb=readSequenceDb):
                print("Divergent base: {}".format(node_to_combine))
                nodes_to_combine.remove(node_to_combine)
                continue

            # We are good to start merging nodes. This requires following steps: (I.) generating new node name,
            # (II.) adding the new node to the subgraph, (III.) adding an edge between the new node and the 2nd degree
            # predecessor, (IV.) adding an edge between the new node and its sucessor and (V.) updating the position of
            # the read in the readPositionDb (VI.) updating the read sequence in readSequenceDb and (VII.) cleanup of
            # old nodes in the subgraph.

            # ad (I.)
            # Generating new header using '|' separating nodes
            # (e.g. '605.2' + '1822.2|979.2|637.1' = '605.2|1822.2|979.2|637.1')
            # Retrieving the value of an overlap (number).
            combined_node = predecessor + '|' + node_to_combine
            overlap_to_predecessor = subgraph[predecessor][node_to_combine]['overlap']

            # ad (II.)
            # Adding a combined node to the networkx graph G
            subgraph.add_node(combined_node)

            # ad (III.)
            # Attching the merged node to its new predecessors (original predecessor's predecessor)
            for predecessors_predecessor in predecessors_predecessors:
                # overlap between predecessors and 2. level predecessors
                overlap = subgraph[predecessors_predecessor][predecessor]['overlap']
                # Adding an edge to the network using combined nodes and overlap value
                subgraph.add_edge(predecessors_predecessor, combined_node, overlap=overlap)

            # ad (IV.)
            # Looping over the successors
            for successor in successors:
                # Retrieving the overlap value betwenn sucessor and the node to combine
                overlap = subgraph[node_to_combine][successor]['overlap']
                # Adding an edge to the network G using combined nodes and the overlap value
                subgraph.add_edge(combined_node, successor, overlap=overlap)

            # ad (V.)
            # Retrieving the read offset, counting from the overlap onwards and adding this to the readPositionDb.
            # Note that reads are only stored separately in the readPositionDb!
            offset = len(readSequenceDb[predecessor]) - overlap_to_predecessor

            for read_id in combined_node.split("|"):
                try:
                    readPositionDb[read_id] += offset
                except:
                    print("readPositionDb: {}".format(readPositionDb))

            # ad (VI.)
            # Updating the readSequenceDb entry

            # readSequenceDb[node_to_combine_new] = readSequenceDb[node_to_combine]

            # Retrieving the sequence of predecessor from the read_database (dictionary of all reads)
            predecessor_seq = readSequenceDb[predecessor]
            # Retrieving the sequence of node to combine from the read database (dictionary of all reads)
            node_seq = readSequenceDb[node_to_combine]
            # Combining the two reads (predecessor + overhang of the sequence node)
            combined_seq = predecessor_seq + node_seq[overlap_to_predecessor:]

            # Adding the combined read to the dictionary of all reads
            readSequenceDb[combined_node] = combined_seq

            ########################################################
            # If we want to clean original nodes fro readSequenceDb
            ########################################################
            if purge_original_nodes == True:
                del readSequenceDb[node_to_combine]
                del readSequenceDb[predecessor]
            ######################################

            # ad (VII.) clean up of old nodes in the subgraph
            # Removing predecessor and node to commbined from the network
            subgraph.remove_node(node_to_combine)
            subgraph.remove_node(predecessor)

            # If node to combine still is in the list of nodes to combine, it's removed
            if node_to_combine in nodes_to_combine:
                nodes_to_combine.remove(node_to_combine)
            # If predecessor is in the list of nodes to combine, it will be removed
            if predecessor in nodes_to_combine:
                nodes_to_combine.remove(predecessor)

    return subgraph, readSequenceDb, readPositionDb


def merge_bifurcation(subgraph, readSequenceDb, readPositionDb, variables):
    # Bifurcation happens in a subgraph due to two parts of sequence being divergent. This needs to be treated
    # delicately and the user must have full control over the merge behaviour.

    if not isclean(subgraph):
        print("Not clean!")

    while True:
        merged = False

        # fork out
        # Starting a collapse_candidate as an empty set (unordered list with unique elements)
        collapse_candidate = set([])

        #try:
        for node in subgraph.nodes():
            # Looping through nodes of the subgraph and attempting to merge nodes that are bifurcating at the very
            # egdge of the subgraph. (I.)_Firstly the program verifies whether the current node is
            # still present in the subgraph as later nodes can be removed by 'merge node' function called below.
            # (II.) Fork is defined as a node with 2 or more successors simultaneously if node doesn't have more then
            # 2 successors, it's skipped (not a bifurcation). (III.) When a bifurcation has been identified, the tip
            # candidates are searched as nodes that have out_degree equal to 0 (no nodes follow up). (IV.) Tip
            # candidates are then removed from the set of bifurcating nodes. If there are subsequently no nodes left,
            # the bifurcation occured at the end of the branch. The tip nodes will be re-examined and the node with
            # the higher number of reads will be passed on. (V.) Nodes with the highest number of
            # reads are selected and passed on the the 'merge_node' function. (VI.) The 'merge_node' function is
            # returning a merged_node and that one is assigned as candidate for a 'graph_collapse' function input.

            # Ad (I.)
            # Searching for nodes that may have been removed in prior loops.
            # Not sure how this ever happens, but it does. Probably for nodes removed during the cycle.
            if node not in subgraph.nodes():
                print("node is (not) in subgraph nodes!")
                continue

            # ad (II.)
            successors = set(subgraph.successors(node))
            # If there is fewer then 2 successors, node is not a bifurcation.
            if len(successors) < 2:
                continue

            # ad (III.)
            # Retrieving the potential tips of the network.
            # Tips are nodes that have no successors.
            # In a network 1-2-3, 1 and 2 have out_degree == 1 and 3 have out_degree == 0

            ############################################
            # This condition is most peculiar. If neither of the successors is a tip candidate, the loop will
            # be terminated. Why would this be the case? Are we only interested in closing branching that is way
            # too short?
            ############################################

            tip_candidates = set([s for s in successors if subgraph.out_degree(s) == 0])
            if len(tip_candidates) == 0:
                # If there are no tip candidates among the branching node's successors, skip the loop.
                print("No tip candidates")
                #import pdb; pdb.set_trace()
                continue

            # ad (IV.)
            # Subtracting sets will remove elements of the subtrahend from minued (yields NOT-successors)
            merge_candidates = successors - tip_candidates

            # ad (V.)
            # Node with the highest number of reads will be searched and passed to the 'merge_node' function.

            # If, after subtraction of tip_candidates there are no dst candidates left, the bifurcation is at
            # the tip only. In the opposite case,

            # A custom function (see below) 'n_reads_in_node' determines the number of reads in the particular node.
            if len(merge_candidates) == 0:
                # Looping through the dst nodes, remove a node with maximum reads
                # n_read_in_node is a custom function splits read header by '|' and counts sequence codes
                toMerge = max([[n_read_in_node(t), t] for t in tip_candidates])[1]
                # Removing the tip candidate from the tip_candidates
                tip_candidates.remove(toMerge)
            else:
                # In this case, we are looping through the internal bifurcation members (excluding tips, which have
                # been removed earlier).
                toMerge = max([[n_read_in_node(d), d] for d in merge_candidates])[1]
                merge_candidates.remove(toMerge)  # remove dst
                # Looping over the merge_candidates list to check for tip nodes (out_degree == 0) after removing
                # the dst_node.
                # TODO: I don't understand how this should ever occur. Check if it ever happens!
                internal_candidates = [d for d in merge_candidates if subgraph.out_degree(d) == 0]
                if len(internal_candidates) > 0:
                    print("Somehow we created a tip here.")
                    print(merge_candidates)
                    #import pdb;pdb.set_trace()
                # Merging the dst_candidates with the tip_candidates
                tip_candidates = tip_candidates.union(merge_candidates)

            if not isclean(subgraph):
                print("Not clean!")

            # Merging the merge_node
            merged_node, subgraph, readSequenceDb, readPositionDb = merge_node(nodesToMergeWith=tip_candidates,
                                                                               nodeToMerge=toMerge,
                                                                               sharedNode=node,
                                                                               subgraph=subgraph,
                                                                               readSequenceDb=readSequenceDb,
                                                                               readPositionDb=readPositionDb,
                                                                               variables=variables,
                                                                               direction=1)

            # If the input isn't empty:
            if merged_node:
                merged = True
                # Adding a node to the collapse candidate set
                collapse_candidate.add(node)

        #except:
        #    print("Error: Checking for nodes in subgraph have failed!")
            #import pdb;
            #pdb.set_trace()
        #    quit()
        if not isclean(subgraph):
            print("Not clean!")
        # ad (VI.)
        # Calling a custom collapse function
        subgraph, readDatabase, readPositionDb = collapse_graph(subgraph=subgraph,
                                                                candidates=list(collapse_candidate),
                                                                readSequenceDb=readSequenceDb,
                                                                readPositionDb=readPositionDb,
                                                                purge_original_nodes=False)

        # fork in
        # Repeating the same procedure from before in the opposite direction.

        collapse_candidate = set([])
        # Looping through nodes in the G network
        for node in subgraph.nodes():
            if not isclean(subgraph):
                print("Not clean!")
        #try:
        for node in subgraph.nodes():

            if not isclean(subgraph):
                print("NOT CLEAN!")
                break

            # ad (I.)
            if node not in subgraph.nodes():
                continue

            # ad (II.) in opposite direction
            predecessors = set(subgraph.predecessors(node))
            if len(predecessors) < 2:
                continue

            # ad (III.)
            tip_candidates = set(
                [p for p in predecessors if subgraph.in_degree(p) == 0])  # and G.out_degree(p) == 1])
            if len(tip_candidates) == 0:
                continue

            # ad (IV.)
            merge_candidates = predecessors - tip_candidates

            if not isclean(subgraph):
                print("Not clean!")

            if not iscleanSet(predecessors):
                print("Not clean!")

            if not iscleanSet(tip_candidates):
                print("Not clean!")

            # ad (V.)
            if len(merge_candidates) == 0:
                # Looping through the dst nodes, remove a node with maximum reads
                # n_read_in_node is a custom function splits read header by '|' and counts sequence codes
                toMerge = max([[n_read_in_node(t), t] for t in tip_candidates])[1]
                # Removing the tip candidate from the tip_candidates
                tip_candidates.remove(toMerge)
            else:
                # Looping through the dst nodes, remove a node with maximum reads
                # Custom function splits read header by '|' and counts sequence codes
                toMerge = max([[n_read_in_node(d), d] for d in merge_candidates])[1]  # only one dst node
                merge_candidates.remove(toMerge)  # remove dst
                # Looping over the dst_candidates list to find if there is a tip_candidate after removing the dst_node
                merge_candidates = [d for d in merge_candidates if subgraph.in_degree(
                    d) == 0]  # and G.out_degree(d) == 1]  # only if its out-deg is 0, a node will be considered tip
                tip_candidates = tip_candidates.union(merge_candidates)

            # Merging the dst_nodes
            merged_node, subgraph, readSequenceDb, readPositionDb = merge_node(nodesToMergeWith=tip_candidates,
                                                                               nodeToMerge=toMerge,
                                                                               sharedNode=node,
                                                                               subgraph=subgraph,
                                                                               variables=variables,
                                                                               readSequenceDb=readSequenceDb,
                                                                               readPositionDb=readPositionDb,
                                                                               direction=-1)

            # If the input isn't empty:
            if merged_node:
                merged = True
                # Adding a node to the collapse candidate set
                collapse_candidate.add(node)
        #except:
        #    print("Error: Checking for nodes in subgraph have failed!")
        #    quit()

        # ad (VI.)
        # Calling a custom collapse function
        subgraph, readDatabase, readPositionDb = collapse_graph(subgraph=subgraph, \
                                                                candidates=list(collapse_candidate), \
                                                                readSequenceDb=readSequenceDb, \
                                                                readPositionDb=readPositionDb, \
                                                                purge_original_nodes=False)

        # If no nodes were merged (no merged_nodes retrieved from 'merge_node'), exit the loop.
        if merged == False:
            break

    # Attempt onece again to merge nodes using the 'collapse_graph'
    subgraph, readSequenceDb, readPositionDb = collapse_graph(subgraph=subgraph, \
                                                              candidates=[], \
                                                              readSequenceDb=readSequenceDb, \
                                                              readPositionDb=readPositionDb, \
                                                              purge_original_nodes=False)

    return subgraph, readSequenceDb, readPositionDb


def remove_bubble(subgraph, readSequenceDb, readPositionDb, variables):
    while True:
        bubble_removed = False
        # Retrieving nodes (e.g ['1107.2|1611.2|2403.1'])
        all_nodes = subgraph.nodes()
        # Starting a set of collapse candidates (unordered list)
        collapse_candidate = set([])

        # Looping through nodes in all_nodes
        for node in all_nodes:
            # In the above it would be '1107.2|1611.2|2403.1'
            # Should the node not be present, skip the loop
            if node not in subgraph.nodes():
                continue

            successors = [s for s in subgraph.successors(node) if
                          subgraph.in_degree(s) == 1 and subgraph.out_degree(s) == 1]
            # If there are no successors (like in the example above), skip the loop
            if len(successors) < 2:
                continue
            # Dictonary of successors?
            d = {}
            print("Removing bubbles!")

            # Looping through the list of successors
            for successor in successors:
                # Level 2 successors (successors of the successor)
                to_node = subgraph.successors(successor)[0]  # successor has only one successor
                # If level 2 successor isn't in a dictionary, it's added
                if to_node not in d:
                    d[to_node] = []
                # Successor is added to a node
                d[to_node].append(successor)

            # If there are more then one successor coming after a node, we have a bubble!
            for to_node in [n for n in d if len(d[n]) > 1]:
                # Merging the bubble
                new_node, subgraph, readSequenceDb, readPositionDb = merge_node(nodesToMergeWith=d[to_node][1:],
                                                                                nodeToMerge=d[to_node][0],
                                                                                sharedNode=node,
                                                                                subgraph=subgraph,
                                                                                readSequenceDb=readSequenceDb,
                                                                                readPositionDb=readPositionDb,
                                                                                variables=variables,
                                                                                direction=1)

                # If there is a result of merging, toggler is changed and the loop continues from the top
                if new_node:
                    bubble_removed = True
                    # Adding new node to the collapse candidates list
                    collapse_candidate.add(new_node)

        subgraph, readSequenceDb, readPositionDb = collapse_graph(subgraph=subgraph, \
                                                                  candidates=list(collapse_candidate), \
                                                                  readSequenceDb=readSequenceDb, \
                                                                  readPositionDb=readPositionDb, \
                                                                  purge_original_nodes=False)

        # If no bubbles were removed the loop is exited
        if not bubble_removed:
            break

    return subgraph, readSequenceDb, readPositionDb


def merge_node(nodesToMergeWith, nodeToMerge, sharedNode, subgraph, readSequenceDb, readPositionDb, variables, direction):

    ################################
    # Renaming variables
    ################################
    # dst -> nodeToMerge
    # dst_seq -> nodeToMergeSeq
    # dst_overlap -> nodeToMergeOverlap
    # dst_remaining -> nodeToMergeRemaining
    ################################
    # src_list -> downstreamNodes -> nodesToMergeWith
    # src -> downstreamNode
    # src_seq -> downstreamNodeSeq
    # src_overlap -> downstreamNodeOverlap -> nodeToMergeWithOverlap
    # src_remaining -> downstreamNodeRemaining -> downstreamNodeOverhangSeq
    # src_list, dst and shared look like this: '1581.2|3294.2;...'
    # shared -> nodeToMergeWith -> sharedNode
    ################################

    # Number of allowed mismatches.
    # TODO : Hardcoded mismatches

    N_MIS_percentage = 3

    # Database of read names that have been changed
    oldReadNames = {}

    ###########################################################
    # Assumingly this part merges a predecessor
    ###########################################################

    # Retrieving a sequence of a node that is supposed to be merged and its overlap with the successor/predecessor
    # depending on the direction specified by the input variable 'direction'.

    # dst is a header of merged sequences (e.g. '2579.2|2295.1|1042.1')
    # retrieving a sequence belonging to the header (dst_seq) (e.g. 'GATTCA...')
    nodeToMergeSeq = readSequenceDb[nodeToMerge]

    # The two nodes that Overlap of the sequences of the two nodes is retrieved.
    if direction == 1:
        shared_ToMergeOverlap = subgraph[sharedNode][nodeToMerge]['overlap']
        shared_ToMergeOverlapSeq = nodeToMergeSeq[:shared_ToMergeOverlap]
        nodeToMergeOverhangSeq = nodeToMergeSeq[shared_ToMergeOverlap:]
    else:
        shared_ToMergeOverlap = subgraph[nodeToMerge][sharedNode]['overlap']
        shared_ToMergeOverlapSeq = nodeToMergeSeq[-shared_ToMergeOverlap:]
        nodeToMergeOverhangSeq = nodeToMergeSeq[:-shared_ToMergeOverlap]

        new_sequence = nodeToMergeOverhangSeq + readSequenceDb[sharedNode]

    # Starting list of nodes to remove and merge
    to_remove = []
    # List of nodes (e.g. tip candidates)
    to_merge = []

    ###########################################################
    # This part merges successors
    ###########################################################

    for nodeToMergeWith in nodesToMergeWith:

        ############################################################

        # TODO : Determine the effect of this decision
        # Function counting reads in header, splitting by '|'
        # Number of reads in the src input (see above what that is)
        # If nodes have very different number of reads, they will not be merged
        if n_read_in_node(nodeToMergeWith) >= 1.2 * n_read_in_node(nodeToMerge):
            continue

        #############################################################

        # if src in oldReadNames:
        #    print("{} have allready been dealt with.".format(src))
        #    continue

        # Retrieving sequences of each of the listed nodes
        try:
            nodeToMergeWithSeq = readSequenceDb[nodeToMergeWith]
        except:
            print("Current src: {}".format(nodeToMergeWith))
            print("This tip has already been dealt with.")
            return None, subgraph, readSequenceDb


        # OVERLAP OF THE TWO SEQUENCES
        # Determining the overlap (integer) of the node with the node to be merged with
        # Sequence that's overhanging from the node to be merged with

        if direction == 1:
            shared_ToMergeWithOverlap = subgraph[sharedNode][nodeToMergeWith]['overlap']
            shared_ToMergeWithOverlapSeq = nodeToMergeWithSeq[shared_ToMergeWithOverlap:]
            nodeToMergeWithOverhangSeq = nodeToMergeSeq[shared_ToMergeWithOverlap:]
        else:
            shared_ToMergeWithOverlap = subgraph[nodeToMergeWith][sharedNode]['overlap']
            shared_ToMergeWithOverlapSeq = nodeToMergeWithSeq[-shared_ToMergeWithOverlap:]
            nodeToMergeWithOverhangSeq = nodeToMergeWithSeq[:-shared_ToMergeWithOverlap]

        overhangsOverlapLength = min(len(nodeToMergeWithOverhangSeq), len(nodeToMergeOverhangSeq))

        # Expressing the number of mismatches depending on the length of overlap and rounding the output to be
        # exactly comparable in the future if statement (floats are not behaving well in an if statement).

        N_MIS = (N_MIS_percentage/100)*overhangsOverlapLength + 1 ############### TODO: JUST FOR NOW
        N_MIS = np.floor(N_MIS)

        # Number of mismatches
        mis = 0
        # Looping over the bases of the overlapping parts of the nodes that are to be merged
        newOverlapSeq = []

        if direction == 1:
            if n_read_in_node(nodeToMerge) > n_read_in_node(nodeToMergeWith):
                dominant = nodeToMergeOverhangSeq
                recessive = nodeToMergeWithOverhangSeq
            else:
                dominant = nodeToMergeWithOverhangSeq
                recessive = nodeToMergeOverhangSeq

            for d,r in zip(dominant, recessive):
                newOverlapSeq.append(d)
                if d != r:
                    mis += 1
                    if mis > N_MIS:
                        break

            # for i in range(overhangsOverlapLength):
            #     # If the downstreamNode overhang and nodeToMerge overhang don't match add a mismatch
            #     if n_read_in_node(nodeToMerge) > n_read_in_node(nodeToMergeWith):
            #         seqOver.append(nodeToMerge[i])
            #     else:
            #         seqOver.append(nodeToMergeWithOverhangSeq[i])
            #     if shared_ToMergeWithOverlapSeq[i] != nodeToMergeOverlapSeq[i]:
            #         mis += 1
            #         print('mis {}'.format(mis))
            #         # If the number of mismatches is larger then defined above, break from the loop
            #         if mis > N_MIS:
            #             break
        else:
            if n_read_in_node(nodeToMerge) > n_read_in_node(nodeToMergeWith):
                dominant = nodeToMergeOverhangSeq
                recessive = nodeToMergeWithOverhangSeq
            else:
                dominant = nodeToMergeWithOverhangSeq
                recessive = nodeToMergeOverhangSeq

            for d,r in zip(reversed(dominant), reversed(recessive)):
                newOverlapSeq.append(d)
                if d != r:
                    mis += 1
                    if mis > N_MIS:
                        break



        # TODO : This needs a revision. No idea what's the point of TIP_SIZE or why this is relevant with mismatches
        # If the number of mismatches is larger then allowed above
        if mis > N_MIS:
            # If number of reads in node is larger then the tip size
            if n_read_in_node(nodeToMergeWith) < variables.TIP_SIZE:
                to_remove.append(nodeToMergeWith)
            continue

        # Offset of the sequence to overlap
        if direction == 1:
            offset = shared_ToMergeOverlap - shared_ToMergeWithOverlap
        else:
            offset = ((len(nodeToMergeSeq) - shared_ToMergeOverlap) -
                      (len(nodeToMergeWithSeq) - shared_ToMergeWithOverlap))

        if offset < -1000:
            print("Offset too high!!")

        for readID in nodeToMergeWith.split("|"):
            readPositionDb[readID] += offset

        # Appending the new sequence to the list to merge
        to_merge.append(nodeToMergeWith)

    # Should there be nothing to remove or merge and finish the function
    if not to_remove + to_merge:
        return None, subgraph, readSequenceDb, readPositionDb

    # Looping throught the list of nodes to remove
    for n in to_remove:
        # Removing node from the network
        subgraph.remove_node(n)

    if to_merge:
        new_node = nodeToMerge + "|" + "|".join(to_merge)
        subgraph = nx.relabel_nodes(subgraph, {nodeToMerge: new_node}, copy=False)
        readSequenceDb[new_node] = readSequenceDb.pop(nodeToMerge)

        for n in to_merge:
            subgraph.remove_node(n)

        # Changing the label (header) in the network
        # subgraph = nx.relabel_nodes(subgraph, {nodeToMerge: new_node_header}, copy=False)
        #
        # # Adding the new node in the final database (THAT HAS ALREADY BEEN DONE!)
        # readSequenceDb[new_node_header] = readSequenceDb.pop(nodeToMerge)

        # Looping through the list of nodes to merge
        for n in to_merge:
            # Removing a node from the network
            try:
                subgraph.remove_node(n)
            except:
                print("The node %s is not in the digraph." % (n,))

        #return None, subgraph, readSequenceDb, readPositionDb
        return new_node, subgraph, readSequenceDb, readPositionDb
    else:
        return nodeToMerge, subgraph, readSequenceDb, readPositionDb


def n_read_in_node(node):
    # Retreaving the number of reads in a given node
    # read_name_list = node.split(";")[0]
    read_name_list = node.split("|")
    return len(read_name_list)


def remove_isolated_node(subgraph, readSequenceDb, variables):
    for node in subgraph.nodes():
        # There are not in and out degrees to the node and the number of reads present is lower then 5
        # (e.g. there are 3 reads in '1107.2|1611.2|2403.1')
        # OR if the read length is shorter then 105% of the input read lenth (-l option)
        if not subgraph.in_degree(node) and not subgraph.out_degree(node) and \
                (n_read_in_node(node) < 5 or len(readSequenceDb[node]) < variables.READ_LEN * 1.05):
            # Remove the node from the network
            subgraph.remove_node(node)

    return subgraph, readSequenceDb


def get_branching_aid(subgraph_orig):
    # The function used in 'get assemblie' function.

    # The reverse is a graph with the same nodes and edges but with the directions of the edges reversed.
    subgraph = subgraph_orig.reverse(copy=True)
    # Dictionary
    d = {}
    starting_nodes = []

    # Looping through nodes of the reversed network to make a dictionary of nodes and locate starting nodes
    for node_str in subgraph.nodes():
        # Node string '3352.2|838.2|2758.1'
        d[node_str] = set([node_str])
        # If there are not in degrees node is a starting node
        if subgraph.in_degree(node_str) == 0:
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
            successors = subgraph.successors(front)
            # Looping through successors
            for successor in successors:
                # Adds the starting node to the successor dictionary entry
                d[successor] = d[successor].union(d[front])
                if successor not in queue:
                    # Adds the successor to the queue (will be analysed later)
                    queue.append(successor)
    return d


def get_all_path(subgraph, future_nodes, cur_path, paths, variables):
    last_node = cur_path[-1]
    successors = subgraph.successors(last_node)

    # ending node, stop recursion.
    if successors == []:
        paths.append(cur_path)
        return paths
    else:
        if len(successors) > 1:
            # Confidence increment is a custom function
            candidate = sorted([[confidence_increment(cur_path, s, future_nodes[s], variables), s] for s in successors])
            next_node = candidate[-1][1]
        else:
            next_node = successors[0]
        return get_all_path(subgraph, future_nodes, cur_path + [next_node], paths, variables)


def confidence_increment(visited_path, next_node, future_nodes, variables):
    # This has something to do with the number of forward and reverse reads present in the same assembly

    # Starting an empty dictionary 'd' and a weighted_num_pair_end variable
    d, weighted_num_pair_end = {}, 0
    # Looping through a visited path
    for idx, node in enumerate(visited_path):
        # To remain explored
        for read_id in node.split(";")[0].split('|'):
            # Read (e.g. 175.1 = read number 175, forward primer)
            # Read number : number of read e.g. 175
            # directionality : refers to either forward or reverse (e.g. 1 = forward primer)
            readNumber, directionality = read_id.split(".")
            d[readNumber] = len(visited_path) - idx - 1

    for node in future_nodes:
        for read_id in node.split(';')[0].split("|"):
            readNumber, directionality = read_id.split(".")
            if readNumber in d:
                weighted_num_pair_end += 1 * (variables.CONFIDENCE_BASE ** d[readNumber])

    return weighted_num_pair_end


def get_contig(path, subgraph, readDatabase):
    # Retrieve contig from the read database
    try:
        contig = readDatabase[path[0]]
    except:
        print("Something is afoot.")

        ####################################################
        # It seems that the "path" is a non-updated version of the headers. Doesn't contain a correct position value.
        ####################################################

    for idx in range(1, len(path)):
        prev, cur = path[idx - 1], path[idx]
        seq = readDatabase[cur]
        overlap = subgraph[prev][cur]["overlap"]
        contig += seq[overlap:]
    return contig


def get_start_end_pos(path, contig):
    # This function is trying to get a position of the read within the template
    # cm : position within a template
    min_cm_st = float("inf")  # Preselected dummy values for the test.
    max_cm_ed = 0

    startPosition = []
    endPosition = []
    for entry in path:
        positions = entry.split(';')[3].split('=')[1].split(',')
        startPosition.append(positions[0])
        endPosition.append(positions[1])

    minStartPosition = min(startPosition)
    maxEndPosition = max(endPosition)

    # for read_id in [r for r in"|".join(path).split("|") if g.read_position_db[r] >= 0 and (g.read_position_db[r] + g.READ_LEN <= len(contig))]:
    #    min_cm_st = min(min_cm_st, g.cm_pos[read_id][0])
    #    max_cm_ed = max(max_cm_ed, g.cm_pos[read_id][1])

    return minStartPosition, maxEndPosition


def get_assemblie(subgraph, readSequenceDb, variables):
    # Branching aid loops through network, finding tips (starting nodes) and appends nodes that follow them to a dictionary entry
    future_nodes = get_branching_aid(subgraph)
    full_genes = []
    scaffold_candidates = []

    starting_nodes = [n for n in subgraph.nodes() if subgraph.in_degree(n) == 0]
    for node in starting_nodes:
        paths = get_all_path(subgraph, future_nodes, [node], [], variables)
        for path in paths:
            contig = get_contig(path, subgraph, readSequenceDb)
            if len(contig) >= variables.FULL_LENGTH:
                if variables.NEED_DEFLANK:
                    # start_pos = min([g.r_pos[r][0] - g.read_position_db[r] for r in path[0].split("|")])
                    # end_pos = max([len(dat.read_db_original[r]) - g.r_pos[r][1] for r in path[-1].split("|")])
                    # deflanked_contig = contig[st_pos : len(contig) - ed_pos]
                    deflanked_contig = 0
                else:
                    deflanked_contig = contig
                full_genes.append([path, deflanked_contig])
            else:
                startPosition, endPosition = get_start_end_pos(path=path, contig=contig)
                if len(contig) > 120:
                    scaffold_candidates.append([path, startPosition, endPosition, contig])

    return full_genes, scaffold_candidates, readSequenceDb


def check_divergent_bases(predecessor, node, graph, readSequenceDb):
    offset = len(readSequenceDb[node]) - graph[predecessor][node]['overlap']

    offset_node = ' ' * offset + readSequenceDb[node]
    readSequenceDb[predecessor]

    for i, j in zip(offset_node, readSequenceDb[predecessor]):
        if i != j and i != ' ':
            return False

    return True
