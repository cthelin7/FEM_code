import numpy as np
import Elements
import Nodes

g = 1.0
h = -1.0

# uxx + f = 0
# u(1) = 0
# -ux(0) = 0

num_elements = 4
nodes_per_element = 3
num_nodes = nodes_per_element*num_elements - (num_elements - 1)
he = 1/num_elements
node_ids = []
for i in range(0, num_nodes):
    Nodes.Node(i+1)
    node_ids.append(i+1)

ng_id = 1
ng = Nodes.nodes_list[ng_id]

# determine active nodes
active_nodes = []
for node in Nodes.nodes_list:
    if node is not ng:
        active_nodes.append(node)

# assign nodes to elements
k = nodes_per_element - 2
m = 0                                   # factor for shifting up initial node number in each element
for i in range(1, num_elements + 1):   # go through last node in each element
    this_element = Elements.Element(i)
    nodes_to_add = []
    for j in range(0, nodes_per_element):                       # go through each node in the element
        nodes_to_add.append(Nodes.nodes_list[(i+m) + j])
    this_element.add_nodes(nodes_to_add)
    m += k

# initialize K matrix
K = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float16)

# create ID matrix
ID = np.zeros((2, len(Nodes.nodes_list)), dtype=np.float16)
glob_eq_id = 1
for i, node_id in enumerate(node_ids):
    ID[0][i] = node_id
    if node_id is ng_id:
        ID[1][i] = 0
    else:
        ID[1][i] = glob_eq_id
        glob_eq_id += 1

# create IEN matrix (map global node ids to element ids and local node ids
IEN = np.zeros((nodes_per_element, num_elements), dtype=int)
for i in range(0, len(Elements.elements_list)):
    for j in range(0, nodes_per_element):
        IEN[j][i] = Elements.elements_list[i].nodes[j].node_id

# create LM matrix which maps to global equation numbers

for e in Elements.elements_list[0]:
    for a in range(0, nodes_per_element):
        for b in range(0, num_elements):
            # calc ke
            ke = Elements.elements_list[e].ke(he)
        # calc fe
        fe = None;

    # assemble into global arrays
    for a in range(0, nodes_per_element):
        # add to K
        P = LM[a][e]
        if P != 0:
            for b in range(0, num_elements):
                Q = LM[b][e]
                if Q != 0:
                    Kpq = Kpq + ke
            # add to F
            Fp = Fp + fe
d = Fp/Kpq
