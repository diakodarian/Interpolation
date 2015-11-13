# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 16:23:28 2015

Interpolation of 3D functions periodic in y and z directons and non-periodic
in x direction evaluated at Chebyshev nodes.

@author: Diako Darian
"""

from numpy import *
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from FourierInterpolation import FourierInterpolation1D
from ChebyshevInterpolation import ChebInterpXdir


def FourierChebyshevInterpolation2D(f,ff):
    N = ff.shape
    n = f.shape
    indices = []
    for i in range(n[1]):
	index = int(i*(float(N[1]-1)/float(n[1]-1)))
	indices += [index]
	ff[:,index] = ChebInterpXdir(n[0]-1, N[0], x, u, c ,f[:,i], ff[:,index])

    for i in range(N[0]):
	ff[i,:] = FourierInterpolation1D(ff[i,indices[:-1]], ff[i,:])
    return ff	

def FourierChebyshevInterpolation3D(f,ff):
    N = ff.shape
    n = f.shape
    indices = []
    for i in range(n[2]):
	index = int(i*(float(N[2]-1)/float(n[2]-1)))
	indices += [index]
	ff[:,:,index] = FourierChebyshevInterpolation2D(f[:,:,i], ff[:,:,index])

    for i in range(N[0]):
	for j in range(N[1]):
	    ff[i, j, :] = FourierInterpolation1D(ff[i,j,indices[:-1]], ff[i,j,:])
    return ff	



if __name__ == "__main__":

    test = "3D"

    if test == "2D":
	
	m = 5; M = 8 # Number of nodes for the coarse and the finer grid
	n = array([2**m, 2**(m-1)])
	N = array([2**M, 2**(M-1)])
	
	points = arange(n[0]+1)

	u = cos(pi*points/n[0])       # Grid points in x-direction
	v = linspace(0,2*pi,n[1]+1)   # Grid points in y-direction
	
	x = linspace(-1.,1.,N[0])
	y = linspace(0,2*pi,N[1]+1)

	c = ones(n[0]-1)
	c = hstack([.5,c,.5])
	c *= (-1)**points

	# Create the mesh
	T = meshgrid(u,v, indexing='ij')             # Coarse mesh
	X = meshgrid(x,y, indexing='ij')             # Finer mesh
	
	fn = sin(T[0])*cos(T[1]) # Coarse f
	f = sin(X[0])*cos(X[1])  # Finer f

	ff = empty((N[0],N[1]))
	
	ff = FourierChebyshevInterpolation2D(fn,ff)
	
	e = ff - f[:,:-1]
	error = linalg.norm(e, inf)
	print error

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X[0][:,:-1], X[1][:,:-1], ff, rstride=1, cstride=1, cmap=cm.coolwarm,
	    linewidth=0, antialiased=False)
	#ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.title('Barycentric Lagrange interpolation')
	plt.show() 
	
    elif test == "3D":
	m = 5; M = 8 # Number of nodes for the coarse and the finer grid, respectively
	n = array([2**m, 2**(m-1), 2**(m-2)])
	N = array([2**M, 2**(M-1), 2**(M-2)])
	L = array([2., 2*pi, 4*pi])
	
	points = arange(n[0]+1)

	u = cos(pi*points/n[0])       # Grid points in x-direction
	v = linspace(0,2*pi,n[1]+1)   # Grid points in y-direction
	w = linspace(0,4*pi,n[2]+1)   # Grid points in z-direction
	
	x = linspace(-1.,1.,N[0])
	y = linspace(0,L[1],N[1]+1)
        z = linspace(0,L[2],N[2]+1)
        
	c = ones(n[0]-1)
	c = hstack([.5,c,.5])
	c *= (-1)**points
	
	T = meshgrid(u,v,w, indexing='ij')
	X = meshgrid(x,y,z, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])*cos(T[2])
	f = sin(X[0])*cos(X[1])*cos(X[2])
	fv = empty((N[0],N[1],N[2]))

	fv = FourierChebyshevInterpolation3D(fn, fv)
	e = fv - f[:,:-1,:-1] 
	error = linalg.norm(e[:,:,-1], inf)
	assert allclose(fv, f[:,:-1,:-1])
	print "error: ", error

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X[0][:,:-1,-1], X[1][:,:-1,-1], fv[:,:,-1], rstride=1, cstride=1, cmap=cm.coolwarm,
			    linewidth=0, antialiased=False)
	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()		