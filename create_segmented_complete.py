import networkx as nx
import matplotlib.pyplot as plt

# Create a complete graph with 5 nodes
G = nx.complete_graph(5)

# Divide each edge into 3 parts and add new nodes
H = nx.MultiGraph()
node_pos = nx.circular_layout(G)  # Original node positions
for u, v in G.edges():
    x1, y1 = node_pos[u]
    x2, y2 = node_pos[v]
    edge_len = ((x2-x1)**2 + (y2-y1)**2)**0.5
    delta_x = (x2-x1)/3
    delta_y = (y2-y1)/3
    new_nodes_pos = [(x1 + delta_x, y1 + delta_y),
                     (x1 + 2*delta_x, y1 + 2*delta_y)]
    new_nodes = [(u, v, i) for i in range(1, 3)]
    H.add_nodes_from(new_nodes)
    for i, pos in zip(new_nodes, new_nodes_pos):
        node_pos[i] = pos  # Update the position of new nodes
    H.add_edge(u, new_nodes[0], label='1')
    H.add_edge(new_nodes[0], new_nodes[1], label='2')
    H.add_edge(new_nodes[1], v, label='3')
H.add_edges_from(G.edges())

print(H.nodes)

# Draw the segmented-complete graph
color_map = []
n=0
while n < 25:
    #if n % 3 == 0:
    if n == 9:
        color_map.append('tab:blue')
    elif n == 2:
        color_map.append('tab:blue')
    elif n == 3:
        color_map.append('tab:blue')
    elif n == 12:
        color_map.append('tab:blue')
    elif n == 6:
        color_map.append('tab:blue')
    else:
        color_map.append('red')
    n+=1
#while n < 25:
    #color_map.append('red')
    #n+=1
nx.draw_networkx_nodes(H, node_pos, node_size=500, node_color=color_map)
nx.draw_networkx_edges(H, node_pos)
#nx.draw_networkx_edge_labels(H, node_pos, edge_labels={(u, v, i): i for u, v, k, i in H.edges(keys=True)})
#nx.draw_networkx_labels(H, node_pos, {i: str(i) for i in range(1, 6)}, font_size=16, font_color='w')
plt.axis('off')
plt.show()
