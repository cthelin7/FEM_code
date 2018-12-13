import numpy as np
import Elements
import Nodes
import funcs
import math
import matplotlib.pyplot as plt
import setup_funcs
import scipy.linalg as scp

# uxx + f = 0
# u(1) = 0      ------ ng1 = last, x = 1     ----- ng1_id = num_nodes
# -ux(0) = 0    ------ nh = first, x = 0    ----- nh_id = 1

# inputs
g = 0.0
h = -0.0

nodes = [1001]
p_val = [1]
rho = 1.0
# E = 1000000
E = 1
mode_choice = "fixed-fixed"
# mode_choice = "free-fixed"

results = []
for p_v in p_val:
    this_p_results = []
    for node_num in nodes:

        n_el = node_num    # num elements
        # n_per_el = 3    # nodes per element

        P = p_v

        n_shape_funcs = P + 1

        # f_choice = "constant"
        f_choice = "quadratic"
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

        he = 1/float(n_el)

        f_coeff = he/6.0

        # setup B (Bernstein basis functions)

        # call binomial coefficient
        # setup_funcs.bernstein()

        # setup C (extraction operator)
        Ce = setup_funcs.extraction_operator_setup(P, n_el)

        # # call extraction operator definition function
        # setup_funcs.extraction_operator_setup(P, n_el)

        # define quadrature rate
        quad_rate = 3

        int_points, w = setup_funcs.gaussian_quadrature(P, quad_rate)

        # compute knot vector
        knot_v = setup_funcs.knot_vector(P, n_el)

        # compute node locations
        x_locations = setup_funcs.greville_abscissae(P, n_el, knot_v)
        num_nodes = len(x_locations)


        ng1_id = num_nodes
        if mode_choice == "fixed-fixed":
            ng2_id = 1
            w_func = funcs.fixed_fixed_wn
        elif mode_choice == "free-fixed":
            w_func = funcs.free_fixed_wn
            ng2_id = 0

        # ng1 = Nodes.nodes_list[ng1_id]

        nh_id = 1
        # nh = Nodes.nodes_list[nh_id]

        node_ids = []
        active_nodes = []
        for i, x_loc in enumerate(x_locations):
            node = Nodes.Node(i + 1)
            node.add_location(x_loc)
            node_ids.append(i + 1)
            if node.node_id != ng1_id and node.node_id != ng2_id:
                active_nodes.append(node)

        ng1 = Nodes.nodes_list[ng1_id]
        if mode_choice == "fixed-fixed":
            ng2 = Nodes.nodes_list[ng2_id]
        elif mode_choice == "free-fixed":
            pass
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

        # # when P = 3, need an extra 0 in Id matrix
        # if P == 3:
        #     ID[-2][1] = 0.0

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
        K = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float64)

        # initialize M
        M = np.zeros((len(active_nodes), len(active_nodes)), dtype=np.float64)

        #initialize F
        F = []
        for node in active_nodes:
            F.append(0.0)

        all_Ne = []
        all_dNe = []
        all_d2Ne = []

        for e in range(0, n_el):
            this_Ne = []
            this_dNe = []
            this_d2Ne = []

            fe = [0.0]*(P + 1)
            ke = []
            me = []
            for a in range(0, P + 1):
                ke_col = []
                me_col = []
                for b in range(0, P + 1):
                    ke_col.append(0.0)
                    me_col.append(0.0)
                ke.append(ke_col)
                me.append(me_col)

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

                # d2B = []
                # for a in range(0, n_shape_funcs):
                #     d2B.append(setup_funcs.d2_bernstein(P, a, int_points[i]))
                # d2N = []
                # for n in range(0, n_shape_funcs):
                #     sum_d2B = 0.0
                #     for j in range(0, n_shape_funcs):
                #         sum_d2B += Ce[e][n][j] * d2B[j]
                #     d2N.append(sum_d2B)
                # this_d2Ne.append(d2N)

                xz = 0.0
                for a in range(0, P + 1):
                    # sum the xa*Na
                    xe = x_locations[IEN[e][a]]
                    xz += xe*N[a]

                xi = x_locations[i] # ?
                # evaluate fe
                fz = f(xz)


                for a in range(0, P + 1):
                    for b in range(0, P + 1):
                        ke[a][b] += E*dN[a]*dN[b]*(2.0/he)*w[i]
                        me[a][b] += N[a]*rho*N[b]*(he/2.0)*w[i]

                    # fe[a] += N[a] * fz * (he / 2.0) * w[i]
                print str(e) + " " + str(i)
            #     print "fe" + str(fe)
            # print "ke"
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
                            M[LM[e][a] - 1][LM[e][b] - 1] += me[a][b]
                    # F[LM[e][a] - 1] += fe[a]
            # for row in K:
            #     print row
        #
        # print "K"
        # for row in K:
        #     print row
        #
        # print "M"
        # for row in M:
        #     print row
        #
        # print "F"
        # for row in F:
        #     print row

        # F_np = np.asarray([F])
        # F_npt = np.ndarray.transpose(F_np)
        # # print "F = " + str(F_npt)
        # K_np = np.asarray(K)
        # # print "K = " + str(K_np)
        #
        # d = np.linalg.solve(K_np, F_npt)
        # # print "d = " + str(d)
        #
        # # add the known value for the right boundary
        # d = np.append(d, 0.0)
        # # print "d"
        # # print d

        eig_values, eig_vectors = scp.eigh(K, M)

        nat_freqs = np.sqrt(eig_values)
        # for eig in eig_values:
        #     nat_freqs.append(math.sqrt(eig))
        # # setup_funcs.scale_vector(eig_vectors)
        # for i in xrange(len(eig_vectors)):
        #     eig_vectors[i] = setup_funcs.scale_vector(eig_vectors[i])

        eigVectorArrayArray = [0]*10
        for columnNum in range(0,10):
            extractedColumn = eig_vectors[:, columnNum]
            if (columnNum == 0):
                extractedColumn = np.interp(extractedColumn, (extractedColumn.min(), extractedColumn.max()), (0, 1))
            else:
                extractedColumn = np.interp(extractedColumn, (extractedColumn.min(), extractedColumn.max()), (-1, 1))
            extractedColumn = np.insert(extractedColumn, 0, 0)
            extractedColumn = np.append(extractedColumn, 0)
            eigVectorArrayArray[columnNum] = extractedColumn


        error = 0.0
        d_error = 0.0

        x_h = []
        y_h = []
        # for e in range(0, n_el):
        #     for i in range(0, len(int_points)):
        #
        #         xz = 0.0
        #         this_Ne = all_Ne[e][i]
        #         this_dNe = all_dNe[e][i]
        #         u_exact = 0.0
        #         for a in range(0, n_shape_funcs):
        #             # sum the xa*Na
        #             xe = x_locations[IEN[e][a]]
        #             xz += xe * this_Ne[a]
        #
        #         dx_dz = he/2.0
        #         dz_dx = 2.0/he
        #
        #         u_exact = u(xz)
        #
        #         uhe = 0.0
        #         duhe = 0.0
        #         for a in range(0, P + 1):
        #             # uhe += d[a] * this_Ne[int_points[i]]
        #             uhe += d[IEN[e][a]] * this_Ne[a]
        #             # duhe += d[a] * this_dNe[int_points[i]] * dz_dx
        #
        #         diff = u_exact - uhe
        #         # d_diff = du(xz) - duhe
        #         error += diff * diff * dx_dz * w[i]
        #         # d_error += d_diff * d_diff * 0.5 * he * w[i]
        #
        #         for this_z in range(-1, 2):
        #             x_h.append(xz)
        #             y_h.append(uhe)

        wn = []
        wn_error = []
        norm_mode_number = []
        for n in xrange(len(eig_vectors)):
            wn.append(w_func(n, 1.0, E, rho))
            wn.append(math.pi*(n + 1))
            wn_error.append(eig_values[n]/wn[-1])
            norm_mode_number.append(float(n + 1)/len(eig_vectors))

        # print "wn:" + str(wn)
        # print "eig:" + str(eig_values)

        sqrt_error = math.sqrt(error)
        sqrt_d_error = math.sqrt(d_error)
        print ""
        # print "error: " + str(error)
        # print "sqrt_error (displacement error): " + str(sqrt_error)

        # print he
        # print "d_error: " + str(d_error)
        # print "sqrt_d_error (derivative error): " + str(sqrt_d_error)

        # get exact solution values
        x = np.arange(0.0, 1.0, 0.01)
        y = u(x)

        this_p_results.append([P, n_el, sqrt_error, he, num_nodes, wn_error, norm_mode_number])
        print P, n_el, sqrt_error, he, num_nodes, wn_error, norm_mode_number

        # plt.plot(x, y, 'r', x_h, y_h, 'g--')
        # plt.title("f=" + f_choice + ", n=" + str(n_el))
        # plt.xlabel("x")
        # plt.ylabel("u(x)")
        # plt.show()

        plt.plot(norm_mode_number, wn_error, 'r')
        plt.title("frequencies")
        plt.xlabel("normed mode number")
        plt.ylabel("normed freq")
        plt.show()

        x_eig = np.arange(0.0, 1.0, 1.0/len(eigVectorArrayArray[0]))
        plt.plot(x_eig[0:-1], eigVectorArrayArray[0], label='mode 1')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[1], label='mode 2')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[2], label='mode 3')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[3], label='mode 4')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[4], label='mode 5')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[5], label='mode 6')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[6], label='mode 7')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[7], label='mode 8')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[8], label='mode 9')
        plt.plot(x_eig[0:-1], eigVectorArrayArray[9], label='mode 10')
        plt.title("mode_shapes")
        plt.xlabel("x")
        plt.ylabel("mode displacement")
        plt.show()
    results.append(this_p_results)
p2_errors = []
p3_errors = []
he_list = []

for a in range(0, len(nodes)):
    p2_errors.append(results[0][a][2])
    # p3_errors.append(results[1][a][2])
    he_list.append(results[0][a][3])

plt.plot(he_list, p2_errors, 'r',  he_list, p3_errors, 'b')
# plt.title("f=" + f_choice + ", n=" + str(n_el))
plt.xlabel("he")
plt.ylabel("error")
plt.yscale('log')
plt.xscale('log')
plt.show()


plt.plot(nodes, p2_errors, 'r')
# plt.title("f=" + f_choice + ", n=" + str(n_el))
plt.xlabel("nodes")
plt.ylabel("error")
plt.yscale('log')
plt.xscale('log')
plt.show()