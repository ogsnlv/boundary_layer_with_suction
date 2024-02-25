import numpy as np
import matplotlib.pyplot as plt

eps = 1e-7
u_inf = 1
u_0 = 0
K = 1000

def calculate_velocity_field(eps,K,c_Q):
    
    #coefficients
    
    A = np.zeros(K)
    B = np.zeros(K)
    C = np.zeros(K)
    F = np.zeros(K)
    alpha = np.zeros(K)
    beta = np.zeros(K)

    #velocity and velocity's antiderivative
    u = np.zeros(K)
    u_derivative = np.zeros(K)
    previous_f = np.zeros(K)
    previous_u = np.zeros(K)

    # domain borders
    a = 0
    b = 10

    #dimensionless coordinate
    eta = np.linspace(a, b, K, endpoint=True)

    #dimensional step 
    h = (b - a)/ float(K)

    #boundary conditions

    B[0] = 1
    C[0] = 0
    F[0] = u_0 / u_inf

    previous_f[0] = c_Q

    A[-1] = 0
    B[-1] = 1
    F[-1] = 1
    alpha[1] = -C[0] / B[0]
    beta[1] = F[0] / B[0]


    error = 1

    while error > eps:

        #find the coefficients -- straight stroke 
        for k in range(K-1):
            A[k] = 4 - h * previous_f[k]
            B[k] = -8
            C[k] = 4 + h * previous_f[k]
            alpha[k + 1] = -C[k] / (A[k] * alpha[k] + B[k])
            beta[k + 1] = (F[k] - A[k] * beta[k]) / (A[k] * alpha[k] + B[k])
            u[-1] = (F[-1] - A[-1] * beta[-1]) / (B[-1] + A[-1] * alpha[-1])

        #finnd the velocity from the known coefficients -- reverse stroke
        for k in range(K - 2, -1, -1):
            u[k] = alpha[k + 1] * u[k + 1] + beta[k + 1]

        error = np.max(np.abs(previous_u - u))

        for k in range(1, K - 1):
            previous_f[k] = np.trapz(u[:k], np.linspace(0, eta[k], k, endpoint=True))

        for k in range(K):
            previous_u[k] = u[k]
        
        u_derivative = np.gradient(u, h)

    return u, u_derivative, eta

u, u_d, eta = calculate_velocity_field(eps,K,0)


plt.figure()
plt.plot(eta,u)
plt.grid()
plt.show()

