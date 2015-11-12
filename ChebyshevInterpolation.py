# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 10:43:17 2015

Interpolation of non-periodic 1D, 2D or 3D functions evaluated at
Chebyshev nodes.

@author: Diako Darian

"""
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt


def ChebInterpXdir(n, N, x, t, c, f, ff):
    numer = zeros(N)
    denom = zeros(N)
    exact = zeros(N)
    for j in range(n+1):
	xdiff = x-t[j]
	temp = c[j]/xdiff
	numer += temp*f[j]
	denom += temp
	exact[abs(xdiff) < 1E-15] = 1
    ff = numer/denom
    jj = where(exact == 1)[0]
    ii = jj*(float(n+1)/float(N))
    ii = ii.astype(int)
    ff[jj] = f[-(ii+1)]
    return ff

def ChebInterpYdir(n,N,ff):
    numer = zeros(N)
    denom = zeros(N)
    exact = zeros(N)
    for j in range(n+1):
	xdiff = x-t[j]
	temp = c[j]/xdiff
	index = int(j*(float(N-1)/float(n)))
	numer += temp*ff[index]
	denom += temp
	exact[abs(xdiff) < 1E-15] = 1

    fff = numer/denom
    jj = where(exact == 1)[0]
    fff[jj] = ff[jj]
    ff[:] = fff
    return ff    
    
def ChebyshevInterpolation1D(f, ff):
    n = len(f)-1
    N = len(fv)
    numer = zeros(N)
    denom = zeros(N)
    exact = zeros(N)
    for j in range(n+1):
	xdiff = x-t[j]
	temp = c[j]/xdiff
	numer += temp*f[j]
	denom += temp
	exact[abs(xdiff) < 1E-15] = 1

    ff[:] = numer/denom
    jj = where(exact == 1)[0]
    ii = jj*(float(n+1)/float(N))
    ii = ii.astype(int)
    ff[jj] = f[-(ii+1)]
    return ff

def ChebyshevInterpolation2D(f, ff):
    n = f.shape[0]-1
    N = ff.shape[0]
    for i in range(n+1):
	index = int(i*(float(N-1)/float(n)))
	ff[:,index] = ChebInterpXdir(n,N,f[:,i],ff[:,index])
    for i in range(N):
	ff[i,:] = ChebInterpYdir(n,N,ff[i,:])
    return ff

def ChebyshevInterpolation3D(f, ff):
    n = f.shape[0]-1
    N = ff.shape[0]
    indices = []
    for i in range(n+1):
	index = int(i*(float(N-1)/float(n)))
	indices += [index]
	fv[:,:,index] = ChebyshevInterpolation2D(f[:,:,i], fv[:,:,index])
    for i in range(N):
	for j in range(N):
	    ff[i,j,:] = ChebInterpYdir(n,N,ff[i,j,:])
    return ff


if __name__ == "__main__":
    
    test = "3D"
    n = 2**5 # Number of nodes for the coarse grid
    N = 2**8 # Number of nodes for the finer grid
    points = arange(n+1)
    t = cos(pi*points/n)     # Coarse grid
    x = linspace(-1., 1., N) # Finer grid
    
    c = ones(n-1)
    c = hstack([.5,c,.5])
    c *= (-1)**points

    if test == "1D":
	fn = sin(t)
	f = sin(x)
	fv = empty(N)
	
	fv = ChebyshevInterpolation1D(fn, fv)
	
	e = fv - f
        error = linalg.norm(e, inf)
	print "error: ", error
	plt.figure()
        plt.plot(x, fv, '.g', t, fn, '.r')
        plt.show()

	
    elif test == "2D":
	T = meshgrid(t, t, indexing='ij')
	X = meshgrid(x, x, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])
	f  = sin(X[0])*cos(X[1])
	fv = zeros((N,N))

	fv = ChebyshevInterpolation2D(fn, fv)
	e = fv - f
	error = linalg.norm(e, inf)
	print "error: ", error
	
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X[0], X[1], fv, rstride=1, cstride=1, cmap=cm.coolwarm,
			    linewidth=0, antialiased=False)
	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()
	
    elif test == "3D":
	T = meshgrid(t,t,t, indexing='ij')
	X = meshgrid(x,x,x, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])*cos(T[2])
	f = sin(X[0])*cos(X[1])*cos(X[2])
	fv = empty((N,N,N))

	fv = ChebyshevInterpolation3D(fn, fv)
	e = fv - f 
	error = linalg.norm(e[:,:,-1], inf)
	assert allclose(fv, f)
	print "error: ", error


	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X[0][:,:,-1], X[1][:,:,-1], fv[:,:,-1], rstride=1, cstride=1, cmap=cm.coolwarm,
			    linewidth=0, antialiased=False)
	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()	