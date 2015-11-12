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
    indices = []
    for i in range(n+1):
	index = int(i*(float(m-1)/float(n)))
	indices += [index]
	ff[:,index] = ChebInterpXdir(n, m, x, u, c ,f[:,i], ff[:,index])

    for i in range(m):
	ff[i,:] = FourierInterpolation1D(ff[i,indices[:-1]], ff[i,:])
    return ff	

def FourierChebyshevInterpolation3D(f,ff):
    indices = []
    for i in range(n):
	index = int(i*(float(m-1)/float(n)))
	indices += [index]
	ff[:,:,index] = FourierChebyshevInterpolation2D(f[:,:,i], ff[:,:,index])

    for i in range(m):
	for j in range(m):
	    ff[i, j, :] = FourierInterpolation1D(ff[i,j,indices], ff[i,j,:])
    return ff	



if __name__ == "__main__":

    test = "3D"
    n = 2**6 # Number of nodes for the coarse grid
    m = 2**9 # Number of nodes for the finer grid

    points = arange(n+1)

    u = cos(pi*points/n)       # Grid points in x-direction
    v = linspace(0,2*pi,n+1)   # Grid points in y-direction
    x = linspace(-1.,1.,m)
    y = linspace(0,2*pi,m+1)
    
    c = ones(n-1)
    c = hstack([.5,c,.5])
    c *= (-1)**points
    
    if test == "2D":
	# Create the mesh
	T = meshgrid(u,v, indexing='ij')             # Coarse mesh
	X = meshgrid(x,y, indexing='ij')             # Finer mesh
	
	fn = sin(T[0])*cos(T[1]) # Coarse f
	f = sin(X[0])*cos(X[1])  # Finer f

	ff = empty((m,m))
	
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
	T = meshgrid(u,v,v, indexing='ij')
	X = meshgrid(x,y,y, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])*cos(T[2])
	f = sin(X[0])*cos(X[1])*cos(X[2])
	fv = empty((m,m,m))

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