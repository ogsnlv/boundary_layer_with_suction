import matplotlib.pyplot as plt
import numpy as np
import os

current_directory = os.path.dirname(os.path.abspath(__file__))




#Case 1: grid size sensitivity

filepath = 1
i = 1
# filepath =  '/Users/levonoganesyan/Desktop/Projects/GitHub/CP1/results/grid_step_size_sensitivity_study/testing_K='+ str(i*1000)+'.txt'

plt.figure()

while True:

    try:
        filepath =  '/Users/levonoganesyan/Desktop/Projects/GitHub/CP1/results/grid_step_size_sensitivity_study/testing_K='+ str(i*1000)+'.txt'

        eta, u= np.loadtxt(filepath).T

        plt.plot(eta,u)
        i+=1
    except:
        FileNotFoundError
        break
    

# plt.plot(N,error)
# plt.xlabel('N')
# plt.ylabel('Error')

# plt.plot(eta,u)
plt.grid()
plt.show()