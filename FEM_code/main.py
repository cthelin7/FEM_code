import numpy as np
import Elements
import Nodes
import funcs
import math
import matplotlib.pyplot as plt
import setup_funcs

# uxx + f = 0
# u(1) = 0      ------ ng = last, x = 1     ----- ng_id = num_nodes
# -ux(0) = 0    ------ nh = first, x = 0    ----- nh_id = 1

# inputs
g = 0.0
h = -0.0

n_el = 3    # num nodes
n_per_el = 3    # nodes per element

P = 2

n_shape_funcs = P + 1

# f_choice = "constant"
f_choice = "linear"
# f_choice = "quadratic"

# define f(x)
if f_choice == "constant":
    f = funcs.constant
    fab = funcs.fe
    u = funcs.exact_constant
    du = funcs.d_exact_constant
    uh_f = funcs.approx_constant
    f_string = ["----", "----"]
elif f_choice == "linear":
    f = funcs.linear
    fab = funcs.fe
    u = funcs.exact_linear
    du = funcs.d_exact_linear
    uh_f = funcs.approx_linear
    f_string = ["x", ""]
elif f_choice == "quadratic":
    f = funcs.quadratic
    fab = funcs.fe
    u = funcs.exact_quadratic
    du = funcs.d_exact_quadratic
    uh_f = funcs.approx_quadratic
    f_string = ["----", "----"]


# num_nodes = n_per_el * n_el - (n_el - 1)
he = 1/float(n_el)

f_coeff = he/6.0

# setup B (Bernstein basis functions)

# call binomial coefficient
# setup_funcs.bernstein()

# setup C (extraction operator)
Ce = setup_funcs.extraction_operator_setup(P, n_el)

# call extraction operator definition function
setup_funcs.extraction_operator_setup(P, n_el)

# define quadrature rate
quad_rate = 3

int_points, int_weights = setup_funcs.gaussian_quadrature(P, quad_rate)


# compute node locations
# compute knot vector
knot_v = setup_funcs.knot_vector(P, n_el)
x_locations = setup_funcs.greville_abscissae(P, n_el, knot_v)
num_nodes = len(x_locations)


ng_id = num_nodes
# ng = Nodes.nodes_list[ng_id]

nh_id = 1
# nh = Nodes.nodes_list[nh_id]

node_ids = []
active_nodes = []
for i, x_loc in enumerate(x_locations):
    node = Nodes.Node(i + 1)
    node.add_location(x_loc)
    node_ids.append(i + 1)
    if node.node_id != ng_id:
        active_nodes.append(node)

ng = Nodes.nodes_list[ng_id]
nh = Nodes.nodes_list[nh_id]


#
#
#
#
# # # add locations
# # for node_id in node_ids:
# #     node = Nodes.nodes_list[node_id]
# #     node.add_location(x_locations)
#
#
# # determine active nodes
# active_nodes = []
# for i, node in enumerate(node_ids):
#     if node != ng_id:
#         active_nodes.append(node)



# create ID matrix
ID = np.zeros((len(Nodes.nodes_list),2), dtype=np.float16).tolist()
glob_eq_id = 1
for n, node_id in enumerate(node_ids):
    ID[n][0] = node_id
    if node_id == ng_id:
        ID[n][1] = 0
    else:
        ID[n][1] = glob_eq_id
        glob_eq_id += 1
# for n, node_id in enumerate(node_ids):
#     ID[0][n] = node_id
#     if Nodes.nodes_list[node_id] == active_nodes[n]:
#         ID[1][n] = glob_eq_id
#         glob_eq_id += 1
#     else:
#         ID[1][n] = 0


# define IEN (map global node ids to element ids and local node ids
IEN = np.zeros((n_el, n_shape_funcs), dtype=np.int).tolist()
for e in range(0, n_el):
    for a in range(0, n_shape_funcs):
        IEN[e][a] = e + a

# define LM (map to global equation numbers)
LM = np.zeros((n_el, n_shape_funcs), dtype=np.int).tolist()
for e in range(0, n_el):
    for a in range(0, n_shape_funcs):
        # for i in range(0, len(ID[0])):
        #     if IEN[a][e] == ID[0][i]:
        #         # if ID[1][i] == num_nodes:
        #         #     print "HEY"
        #         LM[a][e] = ID[1][i]        # likely a better numpy way to do this
        LM[e][a] = ID[IEN[e][a]][1]



# initialize K
K = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float16).tolist()

#initialize F
F = []
for node in active_nodes:
    F.append(0.0)



for e in range(0, n_el):
    # B = []
    # N = [0.0]*n_shape_funcs
    big_N = []
    for i in range(0, len(int_points)):
        B = []
        for a in range(0, n_shape_funcs):
            B.append(setup_funcs.bernstein(P, a, int_points[i]))
        N = []
        for i in range(0, n_shape_funcs):
            sum_B = 0.0
            for j in range(0, n_shape_funcs):
                sum_B += Ce[e][i][j]*B[j]
            N.append(sum_B)
        big_N.append(N)





for e in range(0, n_el):
    ke = Elements.elements_list[e].ke(he).tolist()
    x1 = Nodes.nodes_list[IEN[0][e]].location
    x2 = Nodes.nodes_list[IEN[1][e]].location

    fe = [0, 0]
    fe = fab(f_coeff, f(x1), f(x2))

    # if e == n_el - 1:
    #     print "last one"
    # print e, fe
    # assemble into global arrays
    for a in range(0, n_per_el):
        # add to K
        P = LM[a][e]
        if P != 0:
            for b in range(0, n_per_el):
                Q = LM[b][e]
                if Q != 0:
                    K[P-1][Q-1] = K[P-1][Q-1] + ke[a][b]  # need to loop over ke?

            # # Compute boundary conditions (none in this example)
            # if P == nh_id:
            #     if a == 0:
            #         fe[a] += h
            #     elif a == n_per_el - 1:
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
for e in range(0, n_el):
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
plt.title("f=" + f_choice + ", n=" + str(n_el))
plt.xlabel("x")
plt.ylabel("u(x)")
plt.show()
