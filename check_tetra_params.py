import numpy as np
from cmath import exp
from utility import *
tetrahegrons_in_cartesian = {'tetra0': (np.sqrt(8/9),0,-1/3),
                             'tetra1': (-np.sqrt(2/9), np.sqrt(2/3), -1/3),
                             'tetra2': (-np.sqrt(2/9),-np.sqrt(2/3),-1/3),
                             'tetra3': (0,0,-1)}

tetrahegrons_states = {'tetra0': np.array([[1/np.sqrt(3)], [np.sqrt(2/3)]]),
                             'tetra1': np.array([[1/np.sqrt(3)], [np.sqrt(2/3)*exp(1j*2*np.pi/3)]]),
                             'tetra2': np.array([[1/np.sqrt(3)], [np.sqrt(2/3)*exp(1j*4*np.pi/3)]]),
                             'tetra3': np.array([[0], [1]])}

for key,item in tetrahegrons_states.items():
    rho = state2rho(item)
    bloch = rho2bloch(rho)
    print(key + ' ' + str(bloch))

# from cartesian to polar
for key,item in tetrahegrons_in_cartesian.items():
    polar = cartesian_to_polar(item)
    print(key + ' ' + str(polar))

# from polar to unitary