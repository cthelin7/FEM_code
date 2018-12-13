import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.linalg as sp

elemVector = [1000]
PVector = [1, 2, 3]
eigValueVectorVector = [[0], [0], [0]]
frequencyVectorVector = [[0], [0], [0]]
errorVectorVector = [[0], [0], [0]]
normalizedModeNumberVectorVector = [[0], [0], [0]]
eigVectorArrayArray = [[[0]] * 10, [[0]] * 10, [[0]] * 10]
nodeVectorVector = [[0], [0], [0]]

for PValue in range(0, len(PVector)):
    for elemValue in range(0, len(elemVector)):

        ### Input
        P = PVector[PValue]
        print('P: ' + str(P))
        numElem = elemVector[elemValue] - P + 1
        hElem = 1.0 / numElem
        E = 1
        rho = 1


        ### Setup
        def get_bernstein(P, z):
            bernsteinVector = [None] * (P + 1)
            for a in range(0, P + 1):
                binomCoef = math.factorial(P) / (math.factorial(a) * math.factorial(P - a))
                bernstein = (1 / (2.0 ** P)) * binomCoef * ((1 - z) ** (P - a)) * ((1 + z) ** a)
                bernsteinVector[a] = bernstein
            return bernsteinVector


        def get_bernstein_derivative(P, z):
            bernsteinDerivativeVector = [None] * (P + 1)
            for a in range(0, P + 1):
                binomCoef = math.factorial(P) / (math.factorial(a) * math.factorial(P - a))
                bernsteinDerivative = (1 / (2.0 ** P)) * binomCoef * ((z + 1) ** (a - 1)) * (2 * a + P * (-z - 1)) * (
                (1 - z) ** (-a + P - 1))
                bernsteinDerivativeVector[a] = bernsteinDerivative
            return bernsteinDerivativeVector


        def get_extraction_operator(P, e):
            C = []
            if P == 1:
                C = [[1, 0], [0, 1]]
            if P == 2:
                if numElem == 1:
                    C = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                elif e == 0:
                    C = [[1, 0, 0], [0, 1, 1.0 / 2], [0, 0, 1.0 / 2]]
                elif e == numElem - 1:
                    C = [[1.0 / 2, 0, 0], [1.0 / 2, 1, 0], [0, 0, 1]]
                else:
                    C = [[1.0 / 2, 0, 0], [1.0 / 2, 1, 1.0 / 2], [0, 0, 1.0 / 2]]
            elif P == 3:
                if numElem == 1:
                    C = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
                elif e == 0:
                    C = [[1, 0, 0, 0], [0, 1, 1.0 / 2, 1.0 / 4], [0, 0, 1.0 / 2, 7.0 / 12], [0, 0, 0, 1.0 / 6]]
                elif e == 1:
                    C = [[1.0 / 4, 0, 0, 0], [7.0 / 12, 2.0 / 3, 1.0 / 3, 1.0 / 6],
                         [1.0 / 6, 1.0 / 3, 2.0 / 3, 2.0 / 3], [0, 0, 0, 1.0 / 6]]
                elif e == numElem - 2:
                    C = [[1.0 / 6, 0, 0, 0], [2.0 / 3, 2.0 / 3, 1.0 / 3, 1.0 / 6],
                         [1.0 / 6, 1.0 / 3, 2.0 / 3, 7.0 / 12], [0, 0, 0, 1.0 / 4]]
                elif e == numElem - 1:
                    C = [[1.0 / 6, 0, 0, 0], [7.0 / 12, 1.0 / 2, 0, 0], [1.0 / 4, 1.0 / 2, 1, 0], [0, 0, 0, 1]]
                else:
                    C = [[1.0 / 6, 0, 0, 0], [2.0 / 3, 2.0 / 3, 1.0 / 3, 1.0 / 6], [1.0 / 6, 1.0 / 3, 2.0 / 3, 2.0 / 3],
                         [0, 0, 0, 1.0 / 6]]
            return C


        def get_quadrature_points():
            z = np.array([-1 * math.sqrt(3.0 / 5), 0, math.sqrt(3.0 / 5)])
            w = np.array([5.0 / 9, 8.0 / 9, 5.0 / 9])
            return z, w


        knotVector = []
        if P == 1:
            knotVector = np.zeros(numElem + 3)
            knotVector[0] = 0
            for i in range(1, numElem + 2):
                knotVector[i] = hElem * (i - 1)
            knotVector[numElem + 2] = 1
        elif P == 2:
            knotVector = np.zeros(numElem + 5)
            knotVector[0] = 0
            knotVector[1] = 0
            for i in range(2, numElem + 3):
                knotVector[i] = hElem * (i - 2)
            knotVector[numElem + 3] = 1
            knotVector[numElem + 4] = 1
        elif P == 3:
            knotVector = np.zeros(numElem + 7)
            knotVector[0] = 0
            knotVector[1] = 0
            knotVector[2] = 0
            for i in range(3, numElem + 4):
                knotVector[i] = hElem * (i - 3)
            knotVector[numElem + 4] = 1
            knotVector[numElem + 5] = 1
            knotVector[numElem + 6] = 1

        nodeVector = np.zeros(numElem + P)
        for i in range(0, numElem + P):
            sum = 0
            for j in range(i + 1, i + P + 1):
                sum = sum + knotVector[j]
            nodeVector[i] = 1.0 / P * sum
        nodeVectorVector[PValue] = nodeVector

        skippedIndices = 0
        ID = np.zeros(len(nodeVector))
        for i in range(0, len(nodeVector)):
            if nodeVector[i] == 1:
                ID[i] = 0
                skippedIndices = skippedIndices + 1
            else:
                ID[i] = i + 1 - skippedIndices

        IEN = np.zeros((P + 1, numElem), np.int32)
        for i in range(0, numElem):
            for j in range(0, P + 1):
                IEN[j][i] = i + j

        LM = np.zeros((P + 1, numElem), np.int32)
        for i in range(0, numElem):
            for j in range(0, P + 1):
                LM[j][i] = ID[IEN[j][i]]


        def get_f(x):
            f = x ** 2
            return f


        K = np.zeros((len(nodeVector) - skippedIndices, len(nodeVector) - skippedIndices))
        M = np.zeros((len(nodeVector) - skippedIndices, len(nodeVector) - skippedIndices))
        z, w = get_quadrature_points()
        print('DOF: ' + str(len(K)))

        ### Solve
        for e in range(0, numElem):
            ke = np.zeros((P + 1, P + 1))
            me = np.zeros((P + 1, P + 1))
            for i in range(0, len(z)):
                B = get_bernstein(P, z[i])
                B = np.asarray(B)
                C = get_extraction_operator(P, e)
                C = np.asarray(C)
                N = np.matmul(C, B)

                BPrime = get_bernstein_derivative(P, z[i])
                BPrime = np.asarray(BPrime)
                NPrime = np.matmul(C, BPrime)

                for a in range(0, P + 1):
                    for b in range(0, P + 1):
                        ke[a][b] = ke[a][b] + E * NPrime[a] * NPrime[b] * 2 / hElem * w[i]
                        me[a][b] = me[a][b] + rho * N[a] * N[b] * hElem / 2 * w[i]

            for a in range(0, P + 1):
                if LM[a][e] > 0:
                    for b in range(0, P + 1):
                        if LM[b][e] > 0:
                            K[LM[a][e] - 1][LM[b][e] - 1] = K[LM[a][e] - 1][LM[b][e] - 1] + ke[a][b]
                            M[LM[a][e] - 1][LM[b][e] - 1] = M[LM[a][e] - 1][LM[b][e] - 1] + me[a][b]

        eigValueVector, eigVectorArray = sp.eigh(K, M)
        eigValueVector = np.sqrt(eigValueVector)
        eigValueVectorVector[PValue] = eigValueVector

        for columnNum in range(0, 10):
            extractedColumn = eigVectorArray[:, columnNum]
            if (columnNum == 0):
                extractedColumn = np.interp(extractedColumn, (extractedColumn.min(), extractedColumn.max()), (1, 0))
            else:
                extractedColumn = np.interp(extractedColumn, (extractedColumn.min(), extractedColumn.max()), (-1, 1))
            extractedColumn = np.append(extractedColumn, 0)
            eigVectorArrayArray[PValue][columnNum] = extractedColumn

        frequencyVector = []
        normalizedModeNumberVector = []
        errorVector = []

        for frequencyNum in range(0, 1000):
            frequency = math.pi * (frequencyNum + 1 - 0.5)
            frequencyVector.append(frequency)
            normalizedModeNumber = (frequencyNum + 1) / 1000.0
            normalizedModeNumberVector.append(normalizedModeNumber)
            errorVector.append(eigValueVector[frequencyNum] / frequencyVector[frequencyNum])

        frequencyVectorVector[PValue] = frequencyVector
        normalizedModeNumberVectorVector[PValue] = normalizedModeNumberVector
        errorVectorVector[PValue] = errorVector

plt.figure(1)
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][0], label='Mode = 1')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][1], label='Mode = 2')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][2], label='Mode = 3')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][3], label='Mode = 4')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][4], label='Mode = 5')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][5], label='Mode = 6')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][6], label='Mode = 7')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][7], label='Mode = 8')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][8], label='Mode = 9')
plt.plot(nodeVectorVector[0], eigVectorArrayArray[0][9], label='Mode = 10')
plt.xlabel('Node Location')
plt.ylabel('Displacement Amplitude')
plt.title('P=1 Mode Shapes, Free-Fixed')
plt.legend(loc='upper left')

plt.figure(2)
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][0], label='Mode = 1')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][1], label='Mode = 2')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][2], label='Mode = 3')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][3], label='Mode = 4')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][4], label='Mode = 5')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][5], label='Mode = 6')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][6], label='Mode = 7')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][7], label='Mode = 8')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][8], label='Mode = 9')
plt.plot(nodeVectorVector[1], eigVectorArrayArray[1][9], label='Mode = 10')
plt.xlabel('Node Location')
plt.ylabel('Displacement Amplitude')
plt.title('P=2 Mode Shapes, Free-Fixed')
plt.legend(loc='upper left')

plt.figure(3)
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][0], label='Mode = 1')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][1], label='Mode = 2')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][2], label='Mode = 3')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][3], label='Mode = 4')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][4], label='Mode = 5')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][5], label='Mode = 6')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][6], label='Mode = 7')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][7], label='Mode = 8')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][8], label='Mode = 9')
plt.plot(nodeVectorVector[2], eigVectorArrayArray[2][9], label='Mode = 10')
plt.xlabel('Node Location')
plt.ylabel('Displacement Amplitude')
plt.title('P=3 Mode Shapes, Free-Fixed')
plt.legend(loc='upper left')

plt.figure(4)
plt.plot(normalizedModeNumberVectorVector[0], errorVectorVector[0], 'r--', label='P=1')
plt.plot(normalizedModeNumberVectorVector[1], errorVectorVector[1], 'g--', label='P=2')
plt.plot(normalizedModeNumberVectorVector[2], errorVectorVector[2], 'b--', label='P=3')
plt.xlabel('Normalized Mode Number (n/N)')
plt.ylabel('Frequency Error (wh/w)')
plt.ylim(bottom=0.99, top=1.23)
plt.title('Free-Fixed')
plt.legend(loc='upper left')
plt.show()
