import numpy as np
from src import constants

def calculate_velocity_field(K, c_Q):

    #dimensionless coordinate array
    eta = np.linspace(constants.START, constants.END, K, endpoint = True)
    # step in eta-direction
    h = (constants.END - constants.START) / (K-1)
    #matrix coefficients
    A = np.zeros(K)
    B = np.zeros(K)
    C = np.zeros(K)
    F = np.zeros(K)

    #propagation coefficients
    alpha = np.zeros(K)
    beta = np.zeros(K)

    #velocity
    u = np.zeros(K)
    #velocity's derivative
    u_prime = np.zeros(K)

    f_prev = np.zeros(K)
    u_prev = np.zeros(K)
    difference_squared = np.zeros(K)

    # conditions on the plate
    B[0] = 1
    C[0] = 0
    F[0] = constants.U_0 / constants.U_INF
    f_prev[0] = c_Q

    #conditions in the freestream
    A[-1] = 0
    B[-1] = 1
    F[-1] = 1
    alpha[1] = -C[0] / B[0]
    beta[1] = F[0] / B[0]

    error = float('inf')

    #iteration process till the error between the current iteration velocity and the previous iteration velocity is less than epsilon
    while error > constants.EPS:
        #tridiagonal method
        #straight propagation
        for k in range(K-1):
            A[k] = 4 - h * f_prev[k]
            B[k] = -8
            C[k] = 4 + h * f_prev[k]
            alpha[k + 1] = -C[k] / (A[k] * alpha[k] + B[k])
            beta[k + 1] = (F[k] - A[k] * beta[k]) / (A[k] * alpha[k] + B[k])

        #finding the velocity component near the freestream first
        u[-1] = (F[-1] - A[-1] * beta[-1]) / (B[-1] + A[-1] * alpha[-1])

        #reverse propagation
        for k in range(K-2, -1, -1):
            u[k] = alpha[k + 1] * u[k + 1] + beta[k + 1]



        error = np.sqrt(np.sum((u - u_prev)**2)/K)

        u_prev[:] = u[:]

        f_prev = [np.trapezoid(u_prev[:i+1],eta[:i+1]) for i in range(len(eta))]


    return (eta,f_prev, u)
