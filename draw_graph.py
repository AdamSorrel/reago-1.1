import networkx as nx


nx.draw(subgraph)

pos=nx.spring_layout(subgraph)

#nx.draw_networkx_edge_labels(G=subgraph, pos=pos, font_size=8)

nx.draw_networkx_labels(G=subgraph, pos=pos,font_size=8)

