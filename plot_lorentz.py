from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ctypes import CDLL, POINTER, c_double
__doc__="""Runs Fortran code to solve Lorentz dynamical system and plots 
attractor."""

if __name__ == '__main__':
    # Retrieve data using Fortran code
    lorentz = CDLL('./lorentz_wrapper.so')
    lorentz.c_butterfly()
    
    data = np.empty((279, 4), dtype="double")
    lorentz.c_assimilate_data(data.ctypes.data_as(POINTER(c_double)))
    
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot(data[:,1], data[:,2], data[:,3])
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.grid(True)
    ax1.set_title('Lorentz Attractor')
    
    plt.show()