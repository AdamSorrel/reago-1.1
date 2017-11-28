import networkx as nx
from math import log
import numpy as np

import json

from bokeh.models import HoverTool, ColumnDataSource, WheelZoomTool, LabelSet
from bokeh.plotting import show, figure

def subgraphPlotBokeh(network, readSequenceDb, externalFiles = False, subgraphFile = None, readSequenceDbFile = None):

    # Importing data from files.
    # Network must be saved in a .gml format and the readSequenceDb must be saved in .json format.
    if externalFiles == True:
        network = nx.read_gml(subgraphFile)

        with open(readSequenceDbFile, 'r') as fp:
            readSequenceDb = json.load(fp)

    # Checking if the network nodes are in long verbose format or just a short one. In the former case, nodes will be
    # renamed to shorter name containing only name of the sequences represented.

    flag = False

    for node in network.nodes():
        if 'read_position' in node:
            flag = True
            break

    if flag == True:
        mapping = {}
        for node in network.nodes():
            mapping[node] = node.split(';')[0]

        network = nx.relabel_nodes(network, mapping)

    # Retrieving tip nodes from the network (nodes that have either no predecessor or no successor)

    endTips = []
    beginningTips = []

    for node in network.nodes():
        if len(list(network.successors(node))) == 0:
            endTips.append(node)
        elif len(list(network.predecessors(node))) == 0:
            beginningTips.append(node)


    maxPathLen = 0
    longestPath = []

    for start in beginningTips:
        for end in endTips:
            try:
                paths = list(nx.shortest_simple_paths(network, start, end))

                pathLength = max(map(lambda x: len(x), paths))

                if pathLength > maxPathLen:
                    maxPathLen = pathLength

                    longestPath = []
                    for entry in paths:
                        longestPath.append(entry)
                elif pathLength == maxPathLen:
                    for entry in paths:
                        longestPath.append(entry)

            except:
                # No path between the two nodes has been found.
                continue

    # Finding the longest path between each of the two edges

    iterPath = iter(longestPath[0])

    positionX = 0

    ypositions = list(range(0, 10000, 10))

    previousNode = longestPath[0][0]

    nodePosition = {}

    next(iterPath)

    # Constructing a path for the longest path. This will serve as a backbone of the graph.

    for node in iterPath:
        seqLen = len(readSequenceDb[previousNode])
        overhang = seqLen - network[previousNode][node]['overlap']
        positionX = positionX + overhang

        positionY = ypositions.pop(0)

        nodePosition[node] = (positionX, positionY)

        previousNode = node

    # Saving nodes in a separate database

    nodeDict = {}

    for key in nodePosition:

        nodeDict[key] = nodePosition[key]

    # Looping through all nodes of the subgraph. Removing nodes that have already been used in the backbone.

    nodesToIncorporate = list(network.nodes())

    incorporatedNodes = []

    for node in nodeDict:
        nodesToIncorporate.remove(node)
        incorporatedNodes.append(node)

    # Constructing positions of all remaining nodes (those that have not been incorporated in the backbone).

    # Setting up counter to stop infinite loops
    counter = 0

    while len(nodesToIncorporate) > 0:

        tempDict = {}
        # print(nodePosition)

        for node in nodeDict:
            counter = counter + 1

            # print(len(nodesToIncorporate))

            for successor in network.successors(node):
                if successor not in incorporatedNodes:
                    # print(successor)

                    x, y = nodeDict[node]

                    xNew = x + len(readSequenceDb[successor]) - network[node][successor]['overlap']

                    yNew = ypositions.pop(0)

                    # if len(list(network.successors(node))) > 1:
                    #    yNew = y + 10
                    #    print("node : {}, successors: {}".format(node, len(list(network.successors(node)))))
                    # else:
                    #    yNew = y + 0

                    tempDict[successor] = (xNew, yNew)
                    nodesToIncorporate.remove(successor)
                    incorporatedNodes.append(successor)

            for predecessor in network.predecessors(node):
                if predecessor not in incorporatedNodes:
                    x, y = nodeDict[node]

                    xNew = x - len(readSequenceDb[predecessor]) + network[predecessor][node]['overlap']

                    yNew = ypositions.pop(0)

                    # if len(list(network.predecessors(node))) > 1:
                    #    yNew = y + 10
                    #    print("node : {}, predecessors: {}".format(node, len(list(network.successors(node)))))
                    # else:
                    #    yNew = y + 0

                    tempDict[predecessor] = (xNew, yNew)
                    nodesToIncorporate.remove(predecessor)
                    incorporatedNodes.append(predecessor)

        for entry in tempDict:
            # 3print("Adding {}".format(entry))
            nodeDict[entry] = tempDict[entry]
            # print(len(nodeDict))

        if counter > 10000:
            print("TOO MANY ITERATIONS")
            break

    layout = {}

    for node in nodeDict:
        layout[node] = np.array([nodeDict[node][0], nodeDict[node][1]], dtype='int')


    nodes, nodes_coordinates = zip(*sorted(layout.items()))
    nodes_xs, nodes_ys = list(zip(*nodes_coordinates))

    def get_nodes_specs(_network, _layout, _seqDict):
        d = dict(xleft=[], xright=[], ybottom=[], ytop=[], widths=[], onhover=[], seq=[])

        nodes, nodes_coordinates = zip(*sorted(_layout.items()))
        d['xleft'], d['ybottom'] = list(zip(*nodes_coordinates))

        for i, (xl, yb, node) in enumerate(zip(d['xleft'], d['ybottom'], nodes)):
            d['xright'].append(d['xleft'][i] + len(_seqDict[node]))
            d['ytop'].append(d['ybottom'][i] + 5)

        for node in nodes:
            d['onhover'].append(node)
            d['widths'].append(len(_seqDict[node]))
            d['seq'].append(_seqDict[node])
        return d

    def get_edges_specs(_network, _layout):
        d = dict(xs=[], ys=[], alphas=[], widths=[], onhover=[])
        weights = [d['overlap'] for u, v, d in _network.edges(data=True)]
        max_weight = max(weights)
        calc_alpha = lambda h: 1 + log(h / max_weight) / log(2)
        calc_width = lambda h: 1 + log(10 * h / max_weight) / log(2)

        # example: { ..., ('user47', 'da_bjoerni', {'weight': 3}), ... }
        for u, v, data in _network.edges(data=True):
            d['xs'].append([_layout[u][0], _layout[v][0]])
            d['ys'].append([_layout[u][1], _layout[v][1]])
            d['alphas'].append(calc_alpha(data['overlap']))
            d['widths'].append(calc_width(data['overlap']))
            d['onhover'].append("Overlap: " + str(data['overlap']))
        return d

    nodes_source = ColumnDataSource(get_nodes_specs(_network=network, _layout=layout, _seqDict=readSequenceDb))

    lines_source = ColumnDataSource(get_edges_specs(network, layout))

    hover = HoverTool(tooltips=[('Info', '@onhover')])
    plot = figure(plot_width=1500, plot_height=500,
                  tools=['tap', 'box_zoom', hover, 'reset'])

    # labels = LabelSet(x='xleft', y='ybottom', text='seq', level='glyph', x_offset=0, y_offset=3, source=nodes_source, render_mode='canvas',text_font_size='10pt')

    # r_circles = plot.circle('xs', 'ys', size='widths', color='blue', level = 'overlay', source=nodes_source)

    r_quad = plot.quad(left='xleft', right='xright', top='ytop', bottom='ybottom', color='grey', alpha=0.5,
                       source=nodes_source)

    # p = plot.renderers.append(TextAnnotation(text='seq', offset='2px', glyph='rect', source=nodes_source))




    r_lines = plot.multi_line('xs', 'ys', line_width='widths', alpha='alphas', color='gray', source=lines_source)

    plot.add_tools(WheelZoomTool())
    plot.add_tools()
    # plot.add_layout(labels)

    show(plot)

    return None