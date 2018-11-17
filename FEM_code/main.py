import numpy as np
import Elements
import Nodes
import funcs
import math
import matplotlib.pyplot as plt
import setup_funcs

# uxx + f = 0
# u(1) = 0      ------ ng1 = last, x = 1     ----- ng1_id = num_nodes
# -ux(0) = 0    ------ nh = first, x = 0    ----- nh_id = 1

# inputs
g = 0.0
h = -0.0

nodes = [1, 10, 100, 1000]
p_val = [2, 3]
h_list = [0.1, 0.01, 0.005, 0.002, 0.001]
E = 1000000

results = []
for p_v in p_val:
    this_p_results = []
    for node_num in nodes:

        n_el = node_num    # num elements
        # n_per_el = 3    # nodes per element

        P = p_v

        n_shape_funcs = P + 1

        b_side = 0.005
        h_side = 0.005

        I = (1.0/12.0)*b_side*h_side**3.0

        # f_choice = "constant"
        f_choice = "quadratic"
        # f_choice = "quadratic"

        # define f(x)
        if f_choice == "constant":
            # f = funcs.constant
            fab = funcs.fe
            # u = funcs.exact_constant
            # du = funcs.d_exact_constant
            # uh_f = funcs.approx_constant
            # f_string = ["----", "----"]
        elif f_choice == "linear":
            # f = funcs.linear
            fab = funcs.fe
            # u = funcs.exact_linear
            # du = funcs.d_exact_linear
            # uh_f = funcs.approx_linear
            # f_string = ["x", ""]
        elif f_choice == "quadratic":
            f = funcs.quadratic
            fab = funcs.fe
            u = funcs.exact_quadratic
            # du = funcs.d_exact_quadratic
            # uh_f = funcs.approx_quadratic
            # f_string = ["----", "----"]

        he = 1/float(n_el)

        f_coeff = he/6.0

        # setup B (Bernstein basis functions)

        # call binomial coefficient
        # setup_funcs.bernstein()

        # setup C (extraction operator)
        Ce = setup_funcs.extraction_operator_setup(P, n_el)

        # define quadrature rate
        quad_rate = 3

        int_points, w = setup_funcs.gaussian_quadrature(P, quad_rate)

        # compute knot vector
        knot_v = setup_funcs.knot_vector(P, n_el)

        # compute node locations
        x_locations = setup_funcs.greville_abscissae(P, n_el, knot_v)
        num_nodes = len(x_locations)

        ng1_id = num_nodes
        ng2_id = num_nodes - 1
        # ng1 = Nodes.nodes_list[ng1_id]

        nh_id = 1
        # nh = Nodes.nodes_list[nh_id]

        Nodes.nodes_list = {}

        node_ids = []
        active_nodes = []
        for i, x_loc in enumerate(x_locations):
            node = Nodes.Node(i + 1)
            node.add_location(x_loc)
            node_ids.append(i + 1)
            if node.node_id != ng1_id and node.node_id != ng2_id:
                active_nodes.append(node)

        ng1 = Nodes.nodes_list[ng1_id]
        ng2 = Nodes.nodes_list[ng2_id]  # beam problem: need an extra 0 in ID matrix
        nh = Nodes.nodes_list[nh_id]



        # create ID matrix
        ID = np.zeros((len(Nodes.nodes_list), 2), dtype=np.float16).tolist()
        glob_eq_id = 1
        for n, node_id in enumerate(node_ids):
            ID[n][0] = node_id
            if node_id == ng1_id or node_id == ng2_id:
                ID[n][1] = 0
            else:
                ID[n][1] = glob_eq_id
                glob_eq_id += 1

        # when P = 3, need an extra 0 in ID matrix
        # if P == 3:
        # ID[-2][1] = 0.0     # ?????

        # define IEN (map global node ids to element ids and local node ids)
        IEN = np.zeros((n_el, n_shape_funcs), dtype=np.int).tolist()
        for e in range(0, n_el):
            for a in range(0, n_shape_funcs):
                IEN[e][a] = e + a

        # define LM (map to global equation numbers)
        LM = np.zeros((n_el, n_shape_funcs), dtype=np.int).tolist()
        for e in range(0, n_el):
            for a in range(0, n_shape_funcs):
                LM[e][a] = ID[IEN[e][a]][1]     # for an element's shape function, choose the correct node in ID

        # initialize K
        K = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float16).tolist()



        #initialize F
        F = []
        for i in range(0, len(active_nodes)):
            F.append(0.0)

        print '\n'
        print "ID: " + str(np.shape(ID))
        print "LM: " + str(np.shape(LM))
        print "ActNodes: " + str(len(active_nodes))
        print "K: " + str(np.shape(K))
        print "F: " + str(np.shape(F))


        all_Ne = []
        all_dNe = []
        all_d2Ne = []

        for e in range(0, n_el):
            this_Ne = []
            this_dNe = []
            this_d2Ne = []

            fe = [0.0]*(P + 1)
            ke = []
            for a in range(0, P + 1):
                ke_col = []
                for b in range(0, P + 1):
                    ke_col.append(0.0)
                ke.append(ke_col)

            for i in range(0, len(int_points)):
                B = []
                for a in range(0, n_shape_funcs):
                    B.append(setup_funcs.bernstein(P, a, int_points[i]))
                N = []
                for n in range(0, n_shape_funcs):
                    sum_B = 0.0
                    for j in range(0, n_shape_funcs):
                        sum_B += Ce[e][n][j]*B[j]
                    N.append(sum_B)
                this_Ne.append(N)

                dB = []
                for a in range(0, n_shape_funcs):
                    dB.append(setup_funcs.d_bernstein(P, a, int_points[i]))
                dN = []
                for n in range(0, n_shape_funcs):
                    sum_dB = 0.0
                    for j in range(0, n_shape_funcs):
                        sum_dB += Ce[e][n][j] * dB[j]
                    dN.append(sum_dB)
                this_dNe.append(dN)

                d2B = []
                for a in range(0, n_shape_funcs):
                    d2B.append(setup_funcs.d2_bernstein(P, a, int_points[i]))
                d2N = []
                for n in range(0, n_shape_funcs):
                    sum_d2B = 0.0
                    for j in range(0, n_shape_funcs):
                        sum_d2B += Ce[e][n][j] * d2B[j]
                    d2N.append(sum_d2B)
                this_d2Ne.append(d2N)

                xz = 0.0
                for a in range(0, P + 1):
                    # sum the xa*Na
                    xe = x_locations[IEN[e][a]]
                    xz += xe*N[a]


                # xi = x_locations[i] # ?
                # evaluate fe
                fz = f(xz, h_side)

                dx_dz = he / 2.0
                dz_dx = 2.0 / he

                for a in range(0, P + 1):
                    for b in range(0, P + 1):
                        ke[a][b] += d2N[a]*E*I*d2N[b]*((2.0/he)**(-3))*w[i]
                        # ke[a][b] += (d2N[a]-dN[a]*dx_dz)*(dx_dz**-2.)*E*I*(d2N[b]-dN[b]*dx_dz)*(1/(dx_dz**2.))*w[i]
                        # ke[a][b] += d2N[a] * E * I * d2N[b] * ((2./he)**(-2.)) * w[i]

                    fe[a] += N[a] * fz * (he / 2.0) * w[i]
                # print str(e) + " " + str(i)
                # print fe
            # for row in ke:
            #     print row
            # print " "

            all_Ne.append(this_Ne)
            all_dNe.append(this_dNe)
            all_d2Ne.append(this_d2Ne)
            # for row in Ce[e]:
                # print row

            for a in range(0, P + 1):
                if LM[e][a] > 0.0:
                    for b in range(0, P + 1):
                        if LM[e][b] > 0.0:
                            K[LM[e][a] - 1][LM[e][b] - 1] += ke[a][b]
                    F[LM[e][a] - 1] += fe[a]
            # for row in K:
            #     print row


        F_np = np.asarray([F])
        F_npt = np.ndarray.transpose(F_np)
        # print "F = " + str(F_npt)
        K_np = np.asarray(K)
        # print "K = " + str(K_np)

        d = np.linalg.solve(K_np, F_npt)
        print "d = " + str(d)

        # add the known value for the right boundary
        d = np.append(d, 0.0)
        d = np.append(d, 0.0)       # add two zeros for beam problem
        # print d

        error = 0.0
        d_error = 0.0

        x_h = []
        y_h = []
        # FIND THE ERROR
        for e in range(0, n_el):
            for i in range(0, len(int_points)):

                xz = 0.0
                # this_Ne = all_Ne[e][i]
                this_dNe = all_dNe[e][i]
                B = []
                for a in range(0, n_shape_funcs):
                    B.append(setup_funcs.bernstein(P, a, int_points[i]))
                N = []
                for n in range(0, n_shape_funcs):
                    sum_B = 0.0
                    for j in range(0, n_shape_funcs):
                        sum_B += Ce[e][n][j] * B[j]
                    N.append(sum_B)

                for a in range(0, P + 1):
                    # sum the xa*Na
                    xe = x_locations[IEN[e][a]]
                    # xz += xe * this_Ne[a]
                    xz += xe * N[a]

                dx_dz = he/2.0
                dz_dx = 2.0/he

                u_exact = u(xz, h_side)
                # u_exact = (10*h_side**3.0)/(8*E*I)
                # x_h.append(xz)
                # y_h.append(u_exact)

                uhe = 0.0
                duhe = 0.0
                for a in range(0, P + 1):
                    # uhe += d[a] * this_Ne[int_points[i]]
                    # uhe += d[IEN[e][a]] * this_Ne[a]
                    uhe += d[IEN[e][a]] * N[a]
                    # uhe += d[IEN[e][a]] * this_dNe[a]
                    # duhe += d[a] * this_dNe[int_points[i]] * dz_dx

                diff = u_exact - uhe
                # d_diff = du(xz, h_side) - duhe
                error += diff * diff * dx_dz * w[i]
                # d_error += d_diff * d_diff * 0.5 * he * w[i]

                # for this_z in range(-1, 2):
                x_h.append(xz)
                y_h.append(uhe)



        sqrt_error = math.sqrt(error)
        # sqrt_d_error = math.sqrt(d_error)
        print ""
        # print "error: " + str(error)
        # print "sqrt_error (displacement error): " + str(sqrt_error)

        # print he
        # print "d_error: " + str(d_error)
        # print "sqrt_d_error (derivative error): " + str(sqrt_d_error)

        # get exact solution values
        x = np.arange(0.0, 1.0, 0.01)
        y = u(x, h_side)

        this_p_results.append([P, n_el, sqrt_error, he, num_nodes])
        print P, n_el, sqrt_error, he, num_nodes

        # plt.plot(x, y, 'r', x_h, y_h, 'g--')
        # plt.title("P=" + str(P) + ", n=" + str(n_el))
        # plt.xlabel("x")
        # plt.ylabel("u(x)")
        # plt.show()
    results.append(this_p_results)
p2_errors = []
p3_errors = []
he_list = []

for a in range(0, len(nodes)):
    p2_errors.append(results[0][a][2])
    p3_errors.append(results[1][a][2])
    he_list.append(results[0][a][3])

plt.plot(he_list, p2_errors, 'r',  he_list, p3_errors, 'b')
# plt.title("f=" + f_choice + ", n=" + str(n_el))
plt.xlabel("he")
plt.ylabel("error")
plt.yscale('log')
plt.xscale('log')
plt.show()

#
# plt.plot(nodes, p2_errors, 'r',  nodes, p3_errors, 'b')
# # plt.title("f=" + f_choice + ", n=" + str(n_el))
# plt.xlabel("nodes")
# plt.ylabel("error")
# plt.yscale('log')
# plt.xscale('log')
# plt.show()