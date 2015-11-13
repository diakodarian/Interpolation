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

def ChebInterpYdir(n,N,x,t,c,ff):
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
    n = f.shape
    N = ff.shape
    for i in range(n[1]):
	index = int(i*(float(N[1]-1)/float(n[1]-1)))
	ff[:,index] = ChebInterpXdir(n[0]-1,N[0],xx,x,cx,f[:,i],ff[:,index])
    for i in range(N[0]):
	ff[i,:] = ChebInterpYdir(n[1]-1,N[1],yy,y,cy,ff[i,:])
    return ff

def ChebyshevInterpolation3D(f, ff):
    n = f.shape
    N = ff.shape
    indices = []
    for i in range(n[2]):
	index = int(i*(float(N[2]-1)/float(n[2]-1)))
	indices += [index]
	fv[:,:,index] = ChebyshevInterpolation2D(f[:,:,i], fv[:,:,index])
    for i in range(N[0]):
	for j in range(N[1]):
	    ff[i,j,:] = ChebInterpYdir(n[2]-1,N[2],zz,z,cz,ff[i,j,:])
    return ff


if __name__ == "__main__":
    
    test = "3D"
    
    if test == "1D":
	n = 2**5 # Number of nodes for the coarse grid
        N = 2**8 # Number of nodes for the finer grid
        points = arange(n+1)
        t = cos(pi*points/n)     # Coarse grid
        x = linspace(-1., 1., N) # Finer grid
        
        c = ones(n-1)
        c = hstack([.5,c,.5])
        c *= (-1)**points
    
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
	m = 5; M = 8
	n = array([2**m, 2**(m-1)])
	N = array([2**M, 2**(M-1)])
	
	pointsx = arange(n[0]+1)
	pointsy = arange(n[1]+1)
	
        x = cos(pi*pointsx/n[0])     # Coarse grid
        y = cos(pi*pointsy/n[1])     
        
        cx = ones(n[0]-1)
        cx = hstack([.5,cx,.5])
        cx *= (-1)**pointsx
        
        cy = ones(n[1]-1)
        cy = hstack([.5,cy,.5])
        cy *= (-1)**pointsy
        
        xx = linspace(-1., 1., N[0]) # Finer grid
        yy = linspace(-1., 1., N[1])
        
	T = meshgrid(x,y, indexing='ij')
	X = meshgrid(xx,yy, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])
	f  = sin(X[0])*cos(X[1])
	fv = zeros((N[0],N[1]))

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
	
	m = 5; M = 8
	n = array([2**m, 2**(m-1), 2**(m-2)])
	N = array([2**M, 2**(M-1), 2**(M-2)])
	
	pointsx = arange(n[0]+1)
	pointsy = arange(n[1]+1)
	pointsz = arange(n[2]+1)
	
        x = cos(pi*pointsx/n[0])     # Coarse grid
        y = cos(pi*pointsy/n[1])     
        z = cos(pi*pointsz/n[2])     
        
        cx = ones(n[0]-1)           # Wight function
        cx = hstack([.5,cx,.5])
        cx *= (-1)**pointsx
        
        cy = ones(n[1]-1)
        cy = hstack([.5,cy,.5])
        cy *= (-1)**pointsy
        
        cz = ones(n[2]-1)
        cz = hstack([.5,cz,.5])
        cz *= (-1)**pointsz
        
        xx = linspace(-1., 1., N[0]) # Finer grid
        yy = linspace(-1., 1., N[1]) 
        zz = linspace(-1., 1., N[2]) 
        
	T = meshgrid(x, y, z, indexing='ij')
	X = meshgrid(xx, yy, zz, indexing='ij')
	
	fn = sin(T[0])*cos(T[1])*cos(T[2])
	f = sin(X[0])*cos(X[1])*cos(X[2])
	fv = empty((N[0],N[1],N[2]))

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