import numpy as np
import Elements
import Nodes

g = 0
h = -2

num_elements = 4
num_nodes = num_elements + 1
nodes_per_element = 2

for i in range(0, num_nodes):
    Nodes.Node(i)

ng = Nodes.nodes_list[1]

active_nodes = []
for node in Nodes.nodes_list:
    if node is not ng:
        active_nodes.append(node)

for i in range(nodes_per_element - 1, len(Nodes.nodes_list)):   # go through last node in each element
    this_element = Elements.Element(i-1)
    nodes_to_add = []
    for j in range(0, nodes_per_element):                       # go through each node in the element
        nodes_to_add.append(Nodes.nodes_list[i-j])              # shifts end node numbers up by only 1 for each element
    this_element.add_nodes(nodes_to_add)

K = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float16)


LM = np.zeros((nodes_per_element, num_elements), dtype=int)
for i in range(0,len(Elements.elements_list)):
    for j in range(0, nodes_per_element):
        LM[j][i] = Elements.elements_list[i].nodes[0][j].node_id
# ^^^ not quite right

for e in Elements.elements_list[0]:
    for a in a_list:
        for b in b_list:
            # calc ke
            ke = None;
        # calc fe
        fe = None;

    # assemble into global arrays
    for a in a_list:
        # add to K
        P = LM[a][e]
        if P != 0:
            for b in b_list:
                Q = LM[b][e]
                if Q != 0:
                    Kpq = Kpq + ke
            # add to F
            Fp = Fp + fe
d = Fp/Kpq
