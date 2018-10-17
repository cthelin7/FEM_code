import numpy as np
import Elements
import Nodes
import funcs
import math

g = 0.0
h = -0.0

# uxx + f = 0
# u(1) = 0      ------ ng = last, x = 1     ----- ng_id = num_nodes
# -ux(0) = 0    ------ nh = first, x = 0    ----- nh_id = 1

num_elements = 200
nodes_per_element = 2
num_nodes = nodes_per_element*num_elements - (num_elements - 1)
he = 1/float(num_elements)

# choice = "constant"
choice = "linear"
# choice = "quadratic"

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
for i, node in enumerate(Nodes.nodes_list):
    if node is not ng_id:
        active_nodes.append(node)

if choice == "constant":
    f = funcs.constant
    # fab = funcs.fe_constant
    # f_coeff = he/2.0
    fab = funcs.fe
    u = funcs.exact_constant
    uh_f = funcs.approx_constant
    f_string = ["----", "----"]
elif choice == "linear":
    f = funcs.linear
    # f_coeff = he/6.0
    # fab = funcs.fe_linear
    fab = funcs.fe
    u = funcs.exact_linear
    uh_f = funcs.approx_linear
    f_string = ["x", ""]
elif choice == "quadratic":
    f = funcs.quadratic
    # f_coeff = he/12.0
    # fab = funcs.fe_quadratic
    fab = funcs.fe
    u = funcs.exact_quadratic
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
    # fe[0] = f_coeff * 2.0 * f(x1) + f_coeff * f(x2)
    # fe[1] = f_coeff * f(x1) + f_coeff * 2.0* f(x2)
    fe = fab(f_coeff, f(x1), f(x2))

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

# print F
F_np = np.asarray([F])
# print F_np
F_npt = np.ndarray.transpose(F_np)
print "F = " + str(F_npt)
K_np = np.asarray(K)
# K_np_inv = np.linalg.inv(K)
print "K = " + str(K_np)
# print K_np_inv

# # d = K\F
# d = K_np_inv*F_npt
# print d
d = np.linalg.solve(K_np, F_npt)
print "d = " + str(d)

# const_comp1 = 1.0
# const_comp2 = 0.0
# coeff_comp = 1.0/he
#
# # u_vec = np.zeros((1, num_elements), dtype=np.int).tolist()
# uh = [0, 0]
# for e in range(0, num_elements):
#     const = 0.0
#     coeff = 0.0
#     if LM[0][e] != 0.0:
#         const += const_comp1 * d[e]
#         coeff += -coeff_comp * d[e]
#     if LM[1][e] != 0.0:
#         const += const_comp2*d[e + 1]
#         coeff += coeff_comp*d[e + 1]
#
#     print "uh[" + str(e) + "] = " + str(coeff) + f_string[0] + " + " + str(const) + f_string[1]
#     uh[0] += coeff
#     uh[1] += const
#
#     const_comp1 += 1.0
#     const_comp2 -= 1.0
# print "uh = " + str(uh[0]) + f_string[0] + " + " +  str(uh[1]) + f_string[1]
u_he = []
u_h = [0.0]*num_nodes
# for e in range(0, num_elements):
#     x1 = Nodes.nodes_list[IEN[0][e]].location
#     x2 = Nodes.nodes_list[IEN[1][e]].location
#
#     if LM[0][e] != 0.0:
#         d1 = d[e]
#     else:
#         d1 = 0.0
#     if LM[1][e] != 0.0:
#         d2 = d[e + 1]
#     else:
#         d2 = 0.0
#
#     for a in range(0, nodes_per_element):
#         if a == 0:
#             z = -1.0
#         elif a == 1:
#             z = 1.0
#
#         u_he = ((d1 + d2)/2.0)+(0.5*(x1+x2) + 0.5*z*he)*((d2 - d1)/2.0)
#         u_h[e + a] += u_he
uh = [0, 0]
uh_cs = []
uhz = [0, 0]
uhz_cs = []

for e in range(0, num_elements):
    x1 = Nodes.nodes_list[IEN[0][e]].location
    x2 = Nodes.nodes_list[IEN[1][e]].location

    if LM[0][e] != 0.0:
        d1 = d[e][0]
    else:
        d1 = 0.0
    if LM[1][e] != 0.0:
        d2 = d[e + 1][0]
    else:
        d2 = 0.0

    coeff = ((d2 - d1)/he)
    const = ((d2 - d1)/he)*((-x1 - x2)/2.0) + ((d2 + d1)/2.0)

    print "uh[" + str(e) + "] = " + str(coeff) + f_string[0] + " + " + str(const) + f_string[1]
    uh[0] += coeff
    uh[1] += const
    uh_cs.append([coeff, const])

    coeffz = ((d2 - d1)/2.0)
    constz = ((d2 + d1)/2.0)
    uhz_cs.append([coeffz, constz])

    print "***uhz[" + str(e) + "] = " + str(coeffz) + "z + " + str(constz) + f_string[1]

    uhz[0] += coeffz
    uhz[1] += constz

print uh
print uhz_cs
print uhz


z = [0.0, 0.0, 0.0]
w = [0.0, 0.0, 0.0]

z[0] = -math.sqrt(3.0/5.0)
z[1] = 0.0
z[2] = math.sqrt(3.0/5.0)

w[0] = 5.0/9.0
w[1] = 8.0/9.0
w[2] = 5.0/9.0

print "z = " + str(z)
print "w = " + str(w)

err_sum = 0.0
for e in range(0, num_elements):
    coeffz = uhz_cs[e][0]
    constz = uhz_cs[e][1]

    for i in range(0, 3):
        # sum = (abs(u(z[i]) - uh_f(z[i])))*(abs(u(z[i]) - uh_f(z[i])))*w[i]
        uhz_value = coeffz*z[i] + constz
        convert_to_z = he*z[i] + x1 + x2
        z_cubed = convert_to_z*convert_to_z*convert_to_z
        u_exact = (1.0/6.0)-(1.0/48.0)*z_cubed

        err_temp = abs(u_exact - uhz_value)
        err_temp = err_temp*err_temp*(he/2)*w[i]
        err_sum += err_temp

err = math.sqrt(err_sum)
print err
