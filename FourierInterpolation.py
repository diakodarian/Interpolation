# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 17:03:37 2015

Interpolation of periodic 1D, 2D or 3D functions by using fft.

@author: Diako Darian

"""
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

def evaliptrig(y,N):
    n = len(y)
    if(n%2) == 0:
        c = fft.ifft(y)
        a = zeros(N, dtype=complex )
        a[:n/2] = c[:n/2]
        a[N-n/2:] = c[n/2:]
        v = fft.fft(a);
        return v
    else : raise TypeError , 'odd length'

def FourierInterpolation1D(f, fv):
    N = len(fv)
    fv = real(evaliptrig(f,N))    
    return fv

def FourierInterpolation2D(f, fv):
    N = fv.shape
    n = f.shape 
    indices = []
    for i in range(n[1]):
	index = int(i*(float(N[1]-1)/float(n[1]-1)))
	indices += [index]
	fv[:, index] = real(evaliptrig(f[:, i], N[0]))	
    for i in range(N[0]):
	fv[i, :] = real(evaliptrig(fv[i, indices], N[1]))      
    return fv

def FourierInterpolation3D(f, fv):
    N = fv.shape
    n = f.shape
    indices = []
    for i in range(n[2]):
	index = int(i*(float(N[2]-1)/float(n[2]-1)))
	indices += [index]
	fv[:,:,index] = FourierInterpolation2D(f[:,:,i], fv[:,:,index])	
    for i in range(N[0]):
	for j in range(N[1]):
	    fv[i,j,:] = real(evaliptrig(fv[i,j,indices], N[2]))  
    return fv


if __name__ == "__main__":
    
    test = "1D"
    
    if test == "1D":
	n = 2**5 # Number of nodes for the coarse grid
	N = 2**8 # Number of nodes for the finer grid
	t = linspace(0, 2*pi, n+1) # Coarse grid
	x = linspace(0, 2*pi, N+1) # Finer grid
    
	fn = sin(t)
	f = sin(x)
	fv = empty(N)
	
	fv = FourierInterpolation1D(fn[:-1], fv)
	
	e = fv - f[:-1]
        error = linalg.norm(e, inf)
	print "error: ", error
	plt.figure()
	plt.plot(x[:-1], fv, '.r', x[:-1], f[:-1])
	plt.show()
	
    elif test == "2D":
	m = 5; M = 8
	n = array([2**m, 2**(m-1)])
	N = array([2**M, 2**(M-1)])
	L = array([2*pi, 4*pi])
	
	x = arange(n[0], dtype=float)*L[0]/n[0]
        y = arange(n[1], dtype=float)*L[1]/n[1]

        xx = linspace(0, L[0], N[0]+1) # Finer grid
	yy = linspace(0, L[1], N[1]+1) # Finer grid
	
	T = meshgrid(x, y, indexing='ij')
	X = meshgrid(xx, yy, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])
	ff = sin(X[0])*cos(X[1])
	fv = empty((N[0],N[1]))
	
	fv = FourierInterpolation2D(fn, fv)
	e = fv - ff[:-1,:-1] 
	error = linalg.norm(e, inf)
	print "error: ", error

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X[0][:-1,:-1], X[1][:-1,:-1], fv, rstride=1, cstride=1, cmap=cm.coolwarm,
			    linewidth=0, antialiased=False)
	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()
	
    elif test == "3D":
	m = 5; M = 8
	n = array([2**m, 2**(m-1), 2**(m-2)])
	N = array([2**M, 2**(M-1), 2**(M-2)])
	L = array([2*pi, 4*pi, 2*pi])
	
	x = arange(n[0], dtype=float)*L[0]/n[0]
        y = arange(n[1], dtype=float)*L[1]/n[1]
        z = arange(n[2], dtype=float)*L[2]/n[2]
        
        xx = linspace(0, L[0], N[0]+1) # Finer grid
	yy = linspace(0, L[1], N[1]+1) # Finer grid
	zz = linspace(0, L[2], N[2]+1) # Finer grid	
	
	T = meshgrid(x, y, z, indexing='ij')
	X = meshgrid(xx, yy, zz, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])*cos(T[2])
	ff = sin(X[0])*cos(X[1])*cos(X[2])
	fv = empty((N[0],N[1],N[2]))

	fv = FourierInterpolation3D(fn, fv)
	e = fv - ff[:-1,:-1,:-1] 
	error = linalg.norm(e[:,:,1], inf)
	assert allclose(fv, ff[:-1,:-1,:-1])
	print "error: ", error


	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X[0][:-1,:-1,-1], X[1][:-1,:-1,-1], fv[:,:,-1], rstride=1, cstride=1, cmap=cm.coolwarm,
			    linewidth=0, antialiased=False)
	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()	