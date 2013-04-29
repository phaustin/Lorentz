from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ctypes import CDLL, POINTER, c_double
__doc__="""Runs Fortran code to solve Lorentz dynamical system and plots 
attractor."""

if __name__ == '__main__':
    # Run integrator from Fortran code
    lorentz = CDLL('/Users/Vlad/Dropbox/Fortran/Practice Code/Lorentz/lorentz_wrapper.so')
    data = np.empty((1000, 3), dtype="double")
    lorentz.c_integrate(data.ctypes.data_as(POINTER(c_double)))
    
    # Plot the Lorentz attractor
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot(data[:,0], data[:,1], data[:,2])
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.grid(True)
    ax1.set_title('Lorentz Attractor')
    
    plt.show()