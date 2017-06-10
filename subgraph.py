from celery import Celery

app = Celery('subgraph', broker = 'amqp://localhost', backend='amqp')

@app.subgraph
def subgraph_treatment(subgraph):
    #filename = "subgaph" + str(num) + ".eps"
    #draw_graph(subgraph, filename)
    #num += 1

    full_genes = []
    scaffold_candidates = []
    # Correcting bases that are either 10 times less abundant or under the error_correction_treshold

    correct_sequencing_error(subgraph, 5)

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
