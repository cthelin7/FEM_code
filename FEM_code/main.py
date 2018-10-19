import numpy as np
import Elements
import Nodes
import funcs
import math
import matplotlib.pyplot as plt

g = 0.0
h = -0.0

# uxx + f = 0
# u(1) = 0      ------ ng = last, x = 1     ----- ng_id = num_nodes
# -ux(0) = 0    ------ nh = first, x = 0    ----- nh_id = 1

num_elements = 10000
nodes_per_element = 2
num_nodes = nodes_per_element*num_elements - (num_elements - 1)
he = 1/float(num_elements)

# choice = "constant"
# choice = "linear"
choice = "quadratic"

f_coeff = he/6.0

node_ids = []
for i in range(0, num_nodes):
    node = Nodes.Node(i+1)
    node.add_location(i * he)
    node_ids.append(i+1)


ng_id = num_nodes
ng = Nodes.nodes_list[ng_id]

nh_id = 1
nh = Nodes.nodes_list[nh_id]

# determine active nodes
active_nodes = []
# for i, node in enumerate(Nodes.nodes_list):
for i, node in enumerate(node_ids):
    if node != ng_id:
        active_nodes.append(node)
    # else:
    #     print "GOT THE G NODE"

if choice == "constant":
    f = funcs.constant
    fab = funcs.fe
    u = funcs.exact_constant
    du = funcs.d_exact_constant
    uh_f = funcs.approx_constant
    f_string = ["----", "----"]
elif choice == "linear":
    f = funcs.linear
    fab = funcs.fe
    u = funcs.exact_linear
    du = funcs.d_exact_linear
    uh_f = funcs.approx_linear
    f_string = ["x", ""]
elif choice == "quadratic":
    f = funcs.quadratic
    fab = funcs.fe
    u = funcs.exact_quadratic
    du = funcs.d_exact_quadratic
    uh_f = funcs.approx_quadratic
    f_string = ["----", "----"]

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
    if node_id == ng_id:
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
                # if ID[1][i] == num_nodes:
                #     print "HEY"
                LM[a][e] = ID[1][i]        # likely a better numpy way to do this

for e in range(0, num_elements):
    ke = Elements.elements_list[e].ke(he).tolist()
    x1 = Nodes.nodes_list[IEN[0][e]].location
    x2 = Nodes.nodes_list[IEN[1][e]].location

    fe = [0, 0]
    fe = fab(f_coeff, f(x1), f(x2))

    # if e == num_elements - 1:
    #     print "last one"
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

            # # Compute boundary conditions (none in this example)
            # if P == nh_id:
            #     if a == 0:
            #         fe[a] += h
            #     elif a == nodes_per_element - 1:
            #         fe[a] += -h
            #
            #     if P == ng_id:
            #         if

            # add to F
            F[P-1] = F[P-1] + fe[a]


F_np = np.asarray([F])
F_npt = np.ndarray.transpose(F_np)
# print "F = " + str(F_npt)
K_np = np.asarray(K)
# print "K = " + str(K_np)

d = np.linalg.solve(K_np, F_npt)
# print "d = " + str(d)

# add the known value for the right boundary
d = np.append(d, 0.0)

# print "d = " + str(d)

z = [0.0, 0.0, 0.0]
w = [0.0, 0.0, 0.0]

z[0] = -math.sqrt(3.0/5.0)
z[1] = 0.0
z[2] = math.sqrt(3.0/5.0)

w[0] = 5.0/9.0
w[1] = 8.0/9.0
w[2] = 5.0/9.0

# print "z = " + str(z)
# print "w = " + str(w)

error = 0.0
d_error = 0.0
x_h = []
y_h = []
for e in range(0, num_elements):
    x1 = Nodes.nodes_list[IEN[0][e]].location
    x2 = Nodes.nodes_list[IEN[1][e]].location

    d1 = d[e]
    d2 = d[e + 1]

    dx_dz = he/2.0
    dz_dx = 2.0/he

    for i in range(0, 3):
        xz = x1*funcs.n1_el(z[i]) + x2*funcs.n2_el(z[i])  # the value of x from z values
        uhe = d1*funcs.n1_el(z[i]) + d2*funcs.n2_el(z[i])  # the value of uhe from z values
        duhe = (d1*funcs.dn1_el(z[i]) + d2*funcs.dn2_el(z[i]))*dz_dx

        diff = u(xz) - uhe
        d_diff = du(xz) - duhe
        error += diff*diff*dx_dz*w[i]
        d_error += d_diff*d_diff*0.5*he*w[i]

    for this_z in range(-1, 2):
        x_h.append((he*this_z + x1 + x2)/2.0)
        y_h.append(d1*funcs.n1_el(this_z) + d2*funcs.n2_el(this_z))

sqrt_error = math.sqrt(error)
sqrt_d_error = math.sqrt(d_error)
print ""
print "error: " + str(error)
print "sqrt_error (displacement error): " + str(sqrt_error)

print he
print "d_error: " + str(d_error)
print "sqrt_d_error (derivative error): " + str(sqrt_d_error)

# get exact solution values
x = np.arange(0.0, 1.0, 0.01)
y = u(x)

# get node location values
# x_node_loc = []
# y_node_loc = []
# for node in node_ids:
#     x_node_loc.append(Nodes.nodes_list[node].location)
#     y_node_loc.append(0)

# print x
# print y
# plt.plot(x, y, 'r--', x_h, y_h, 'g--', x_node_loc, y_node_loc, 'bo')
plt.plot(x, y, 'r--', x_h, y_h, 'g--')
plt.title("f="+choice+", n="+str(num_elements))
plt.xlabel("x")
plt.ylabel("u(x)")
plt.show()
