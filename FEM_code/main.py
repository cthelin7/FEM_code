import numpy as np
import Elements
import Nodes
import funcs

g = 1.0
h = -1.0

# uxx + f = 0
# u(1) = 0
# -ux(0) = 0

num_elements = 4
nodes_per_element = 2
num_nodes = nodes_per_element*num_elements - (num_elements - 1)
he = 1/float(num_elements)

# choice = "constant"
choice = "linear"
# choice = "quadratic"

f_coeff = 1.0

if choice == "constant":
    f = funcs.constant
    f_coeff = he/2.0
elif choice == "linear":
    f = funcs.linear
    f_coeff = he/6.0
elif choice == "quadratic":
    f = funcs.quadratic
    f_coeff = he/12.0

node_ids = []
for i in range(0, num_nodes):
    node = Nodes.Node(i+1)
    node.add_location(i * he)
    node_ids.append(i+1)


ng_id = 5
ng = Nodes.nodes_list[ng_id]

# determine active nodes
active_nodes = []
# for node in Nodes.nodes_list:
#     if node is not ng:
#         active_nodes.append(node)

for i, node in enumerate(Nodes.nodes_list):
    if node is not ng_id:
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
K = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float16).tolist()
# F = np.zeros((1, len(active_nodes)), dtype=np.float16)
F = []
for node in active_nodes:
    F.append(0.0)

# create ID matrix
ID = np.zeros((2, len(Nodes.nodes_list)), dtype=np.float16).tolist()
glob_eq_id = 1
for n, node_id in enumerate(node_ids):
    ID[0][n] = node_id
    if node_id is ng_id:
        ID[1][n] = 0
    else:
        ID[1][n] = glob_eq_id
        glob_eq_id += 1

# create IEN matrix (map global node ids to element ids and local node ids
IEN = np.zeros((nodes_per_element, num_elements), dtype=np.int).tolist()
for e in range(0, len(Elements.elements_list)):
    for a in range(0, nodes_per_element):
        IEN[a][e] = Elements.elements_list[e].nodes[a].node_id

# create LM matrix which maps to global equation numbers
LM = np.zeros((nodes_per_element, num_elements), dtype=np.int).tolist()
for e in range(0, len(Elements.elements_list)):
    for a in range(0, nodes_per_element):
        for i in range(0, len(ID[0])):
            if IEN[a][e] == ID[0][i]:
                LM[a][e] = ID[1][i]        # likely a better numpy way to do this

for e in range(0, num_elements):
    # for a in range(0, nodes_per_element):
        # for b in range(0, num_elements):
        #     # calc ke
        #     ke = Elements.elements_list[e].ke(he)
        # calc fe

    ke = Elements.elements_list[e].ke(he).tolist()
    x1 = Nodes.nodes_list[IEN[0][e]].location
    x2 = Nodes.nodes_list[IEN[1][e]].location

    fe = [0, 0]
    fe[0] = f_coeff * 2.0 * f(x1) + f_coeff * f(x2)
    fe[1] = f_coeff * f(x1) + f_coeff * 2.0* f(x2)

    # print e, fe
    # assemble into global arrays
    for a in range(0, nodes_per_element):
        # add to K
        P = LM[a][e]
        if P != 0:
            for b in range(0, nodes_per_element):
                Q = LM[b][e]
                if Q != 0:
                    K[P-1][Q-1] = K[P-1][Q-1] + ke[a][b]  # need to loop over ke?
            # add to F
            F[P-1] = F[P-1] + fe[a]
    # print F
print F
# d = F/K
