# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 17:40:58 2021

@author: llucv
"""

import numpy as np
import matplotlib.pyplot as plt

def norm(z):
    norm=np.sqrt(np.real(z*np.conj(z)))
    return norm


def mandelbrot(n_iter,n_p,x,y,delta):
    c=np.zeros((n_p,n_p,2),dtype=complex)
    
    x_min=x
    x_max=x+delta
    y_min=y
    y_max=y+delta
    x_l=np.linspace(x_min,x_max,n_p)
    y_l=np.linspace(y_min,y_max,n_p)
    
    c_num=x_l[:,np.newaxis]+y_l[np.newaxis,:]*1j
    c_iter=np.zeros((n_p,n_p),dtype=complex)
    
    z=np.zeros((n_p,n_p),dtype=complex)
    zero_mat=np.zeros((n_p,n_p),dtype=complex)
    
    for i in range(n_iter):
        n=i+1
        norm_z=norm(z)
        z=np.where(norm_z>100,z,z*z+c_num)
        norm_z_post=norm(z)
        c_iter=np.where((norm_z>2)|(norm_z_post<=2),
                        zero_mat,n)+c_iter
        print(n)

    
    c_iter=np.where((norm_z>2)|(norm_z_post>2),c_iter,n_iter+1)
    c[:,:,0]=c_num[:,:]
    c[:,:,1]=c_iter[:,:]
    return c


#prova____________________________________________________________
"""
A=np.array([[1,3],[2,1]])
C=np.array([[1,4],[3,3]])
D=np.array([[1,1],[1,1]])
B=np.where((A>2)|(C<=2),D,0)
print(B)

B=np.where((A>2),D,0)
print(B)

B=np.where((C<=2),D,0)
print(B)

A=np.linspace(0,10,4)
B=np.linspace(0,10,4)
C=A[:,np.newaxis]+B[np.newaxis,:]*1j
print(C)

A=np.array([[4,9],[5,4]])
print(A[1,0])
plt.imshow(A.transpose(),origin="lower")
"""
#____________________________________________________________


n_iter=1000
n_p=7500
x=-1.2551
y=-0.3834
delta=0.0002
iters=np.zeros((n_p,n_p))

           #cmap="nipy_spectral,"
Mand=(mandelbrot(n_iter, n_p, x, y, delta))
iters[:,:]=-np.real(Mand[:,:,1])
plt.imshow(iters.transpose(),origin='lower',extent=(x,x+delta,y,y+delta),
           cmap="nipy_spectral",interpolation='gaussian')
plt.axis('off')
plt.savefig(str(x)+'x_'+str(y)+'y_'+str(delta)+'h_'+str(n_p)+'p_'
            +str(n_iter)+'iter_mandelbrot.png',dpi=3500,bbox_inches='tight')



    
    


