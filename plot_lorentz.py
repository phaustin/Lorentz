'''Runs Fortran code to solve Lorentz dynamical system and plots attractor.'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import lorentz

def assimilate_data(file_name):
    """Assimilate initial data from specified file."""
    
    return np.genfromtxt(file_name)

if __name__ == '__main__':
    
    data = assimilate_data('data.txt')
    print data
    
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot(data[:,1], data[:,2], data[:,3])
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.grid(True)
    ax1.set_title('Lorentz Attractor')
    
    # Try and run some Fortran code here
    lorentz.hello()
    
    plt.show()
    