import numpy as np
import os
import unittest
import data
from src import constants, solution
from matplotlib import pyplot as plt


class TestMeshSensitivity(unittest.case.TestCase):
    def test_mesh_sensitivity(self):

        plt.figure()
        # plt.yscale('log')
        plt.grid()

        error = float('inf')
        u_0_prev = float('inf')
        result_K = []
        result_error = []
        K = 0

        #load experimental data
        eta_exp = np.genfromtxt('data/eta_exp.dat').T
        u_exp = np.genfromtxt('data/u_exp.dat').T


        while error > constants.EPS:

            K+=20000

            eta,_,u = solution.calculate_velocity_field(K,0)

            error = abs(u[0] - u_exp[0])


            print(K, error)

            result_K.append(K)
            result_error.append(error)

            if K%10000 == 0:
                plt.plot(eta,u, label=str(K))



        plt.plot(eta_exp, u_exp, label = 'Exp')
        plt.scatter(eta_exp, u_exp, label = 'Exp')
        plt.legend()
        plt.show()

        return K
