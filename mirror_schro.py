# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 11:30:39 2021

@author: llucv
"""

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import matplotlib.animation as ani

hbar=1
pi=np.pi

def psi_0(Nx,x0,p0,dev,dx):
    psi=np.zeros((Nx),dtype=np.complex128)
    x=np.linspace(0,Nx,Nx,endpoint=False)
    
    x=x*dx-dx*(Nx-1)/2
    
    x0=x0*dx-dx*(Nx-1)/2
    
    dev=dev*dx
    
    const=(1/(2*pi*dev**2))**(1/4)
    psi[:]=const*np.exp(-(((x[:]-x0)**2)/(dev**2))/4+\
                               1j*p0*(x[:]-x0/2)/hbar)

    return psi

def pot_wall(Nx,x0,pot_max,dev,dx):
    x=np.linspace(0,Nx,Nx,endpoint=False,dtype=np.complex128)
    x=x*dx-dx*(Nx-1)/2
    dev=dev*dx
    
    V_x=np.zeros((Nx),dtype=np.complex128)
    
    if dev>0:
        x0=x0*dx-dx*(Nx-1)/2
        print(x0)
        V_x[:]=pot_max*np.exp(-(((x[:]-x0)/dev)**2)/2)
    else:
        V_x[x0]=pot_max

    return V_x


@jit(nopython=True)
def tridiag(A,B,C,D):
    N=len(D)
    x=np.zeros((N),dtype=np.complex128)
    
    for i in range(1,N):
        W=A[i]/B[i-1]
        B[i] = B[i] - W * C[i - 1]
        D[i] = D[i] - W * D[i - 1]
    
    x[N-1]=D[N-1]/B[N-1]
    
    for j in range(1,N):
        i=N-1-j
        x[i] = (D[i] - C[i] * x[i + 1]) / B[i]
        
    return x

@jit
def psi_ev_ck(Nx,Nt,psi_0,V,dx,dt):
    r=dt/(4*dx*dx)
    psi_t=np.zeros((Nx,Nt),dtype=np.complex128)
    psi_t[:,0]=psi_0[:]
    
    
    
    A=np.zeros((Nx),dtype=np.complex128)
    B=np.zeros((Nx),dtype=np.complex128)
    C=np.zeros((Nx),dtype=np.complex128)
    D=np.zeros((Nx),dtype=np.complex128)
    
    
    for i in range(Nt-1):
        D[1:-1]=(1-1j*dt*V[1:-1]/2-2*r*1j)*psi_t[1:-1,i]+\
            1j*r*(psi_t[:-2,i]+psi_t[2:,i])
            
        D[0]=(1-1j*dt*V[0]/2-2*r*1j)*psi_t[0,i]+1j*r*(psi_t[1,i])
        
        D[Nx-1]=(1-1j*dt*V[Nx-1]/2-2*r*1j)*psi_t[Nx-1,i]+\
            1j*r*(psi_t[Nx-2,i])
            
        A[1:]=-1j*r
        A[0]=0.
        B[:]=1+1j*dt*V[:]/2+2*r*1j
        C[0:-1]=-1j*r
        C[Nx-1]=0.
        
        psi_t[:,i+1]=tridiag(A,B,C,D)
        
    
    return psi_t


@jit
def prob_dens(psi):
    prob=np.real(psi*np.conj(psi))
    return prob

def trapezis_1D(f,h):
    f_shape=np.shape(f)
    suma=0
    for i in range(f_shape[0]-1):
        suma=(f[i+1]+f[i])*h/2+\
                 suma
        
    return suma

def psi_momentum(Nx,Nt,psi,dx,dt,k_lim):
    x=np.linspace(0,Nx,Nx,endpoint=False,dtype=np.complex128)
    x=x*dx-dx*(Nx-1)/2
    
    k=np.linspace(-k_lim,k_lim,Nx,endpoint=False,dtype=np.complex128)
    
    psi_k=np.zeros((Nx,Nt),dtype=np.complex128)
    integ=np.zeros((Nx,Nt),dtype=np.complex128)
    
    print('u here')
    
    for t in range(Nt):
        for i in range(Nx):
            integ[:,t]=np.exp(-1j*k[i]*x[:])*psi[:,t]/(np.sqrt(2*pi))
            psi_k[i,t]=trapezis_1D(integ[:,t],dx)
            
        print(t)   
            
    return psi_k
            
def expected_values(psi,psi_k,Nx,dx,k_lim):
    prob=prob_dens(psi)
    prob_k=prob_dens(psi_k)
    dk=(2*k_lim)/Nx
    x_value=0
    k_value=0
    
    for i in range(Nx):
        x_value=(dx*i-dx*(Nx-1)/2)*prob[i]*dx+x_value
        k_value=(dk*i-dk*(Nx-1)/2)*prob_k[i]*dk+k_value
        
        
    mean_val=np.zeros(2)
    mean_val[0]=x_value
    mean_val[1]=k_value
    
    return mean_val
    
    
            
    
        
    
Nt=2000
Nx=2001

x0=int(Nx/4)
p0=10
k_lim=10*p0
dk=2*k_lim/Nx
dev=50

dt=1/60
dx=2*dt

x0_barr=int((Nx-1)/2)
dev_barr=1
pot_max=106.79467777726097

prob_x=np.zeros((Nx,Nt))
prob_x_t=np.zeros((Nx))
x=np.linspace(0,Nx,Nx,endpoint=False)
x=x*dx-dx*(Nx-1)/2
psi_0=psi_0(Nx, x0, p0, dev, dx)
V=pot_wall(Nx, x0_barr, pot_max, dev_barr, dx)

psi=psi_ev_ck(Nx, Nt, psi_0, V, dx, dt)
print('evo done')

for t in range(Nt):
    prob_x[:,t]=prob_dens(psi[:,t])
print('prob calc done')


'''    
prob_x_1=np.zeros((Nx))
prob_x_2=np.zeros((Nx))
pot_max_1=50
pot_max_2=1000
for j in range(200):
    psi_0_j=psi_0(Nx, x0, p0, dev, dx)
    
    V1=pot_wall(Nx, x0_barr, pot_max_1, dev_barr, dx)
    V2=pot_wall(Nx, x0_barr, pot_max_2, dev_barr, dx)
    
    psi_1=psi_ev_ck(Nx, Nt, psi_0_j, V1, dx, dt)
    psi_2=psi_ev_ck(Nx, Nt, psi_0_j, V2, dx, dt)
    
    prob_x_1=prob_dens(psi_1[:,Nt-2])
    prob_x_2=prob_dens(psi_2[:,Nt-2])
    
    dif1=-np.max(prob_x_1[:int((Nx-1)/2)])+np.max(prob_x_1[int((Nx-1)/2):])
    print(dif1)
    dif2=-np.max(prob_x_2[:int((Nx-1)/2)])+np.max(prob_x_2[int((Nx-1)/2):])
    print(dif2)
    
    half_pot=(pot_max_2-pot_max_1)/2
    
    if dif2<0:
        pot_max_2=pot_max_1+half_pot
        pot_max_1=pot_max_1
    
    else:
        pot_max_1=pot_max_1+half_pot
        pot_max_2=pot_max_1+2*half_pot
    
    print(pot_max_1)
    print(pot_max_2)
'''
    

max_y=1.5*np.max(prob_x)

'''
def update(frame):
    t=frame
    prob_x_t[:]=prob_x[:,t]
    plt.cla()
    plt.plot(x,prob_x_t)
    plt.plot(x,np.real(V)/20)

    plt.axis([x[0],x[Nx-1],0,max_y])

fig = plt.figure(figsize=(8, 8))
ax1 = plt.subplot()

Writer = ani.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)

anim = ani.FuncAnimation(fig, update, 
                               frames = Nt, 
                               blit = False, interval=1/60)

anim.save('prob_tunnel_effect.mp4', writer=writer)
'''


'''
y_min=np.min(np.real(psi))
y_max=np.max(np.real(psi))
print(y_min,y_max)

z_min=np.min(np.imag(psi))
z_max=np.max(np.imag(psi))
print(z_min,z_max)

fig = plt.figure()

fig.set_size_inches(10, 10)

ax = fig.add_subplot(111,projection='3d')

def update(frame):
    t=int(frame/2)
    print(t)
    
    ax.cla()
    ax.set_xlim(x[0],x[Nx-1])
    ax.set_ylim(y_min,y_max)
    ax.set_zlim(z_min,z_max)
    ax.plot3D(x, np.real(psi[:,t]), np.imag(psi[:,t]))
    ax.scatter3D(0.,0.,0.)
    


Writer = ani.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)

anim = ani.FuncAnimation(fig, update, 
                               frames = 2*Nt, 
                               blit = False, interval=1/60)

anim.save('wave_tunnel_effect.mp4', writer=writer)
'''

Nt_r=int(Nt/3)#Nt/3

A=6

psi_k=psi_momentum(Nx, Nt_r, psi, dx, dt, k_lim)

prob_k=prob_dens(psi_k)
prob_x=prob_dens(psi)

print('probs done')

phase_space_cont=np.zeros((Nx,Nx,2))
phase_space_evo=np.zeros((Nx,Nx))
phase_space_evo[:,:]=np.tensordot(prob_x[:,0],prob_k[:,0],axes=0)
phase_space_cont[:,:,0]=phase_space_evo[:,:]

print('fig1')

fig1 = plt.figure()

ax1 = fig1.add_subplot()

def update1(frame):
    t=int(frame)
    print(t)
    
    for i in range(A):
        phase_space_evo[:,:]=np.tensordot(prob_x[:,t*A+i],
                                          prob_k[:,t*A+i],axes=0)
    
    ax1.imshow(phase_space_evo[:,:].transpose(),origin='lower',
               extent=(-dx*(Nx-1)/2,dx*(Nx-1)/2,-k_lim/10,+k_lim/10)
           ,interpolation='gaussian',aspect='auto')
    


fig1.set_size_inches(8, 8)
    

Writer = ani.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)

anim = ani.FuncAnimation(fig1, update1, 
                               frames = int(Nt_r/A), 
                               blit = False, interval=1/60)

anim.save('phase_space_evo.mp4', writer=writer, dpi=100)

fig2 = plt.figure()

ax2 = fig2.add_subplot()

print('fig2')

def update2(frame):
    t=int(frame)
    print(t)
    
    ax2.imshow(phase_space_cont[:,:,0].transpose(),origin='lower',
               extent=(-dx*(Nx-1)/2,dx*(Nx-1)/2,-k_lim/10,+k_lim/10)
           ,interpolation='gaussian',aspect='auto')
    
    for i in range(A):
        phase_space_cont[:,:,1]=np.tensordot(prob_x[:,t*A+i],
                                             prob_k[:,t*A+i],axes=0)\
                            +phase_space_cont[:,:,0]
        
        phase_space_cont[:,:,0]=phase_space_cont[:,:,1]
    
    

    
    
    
    
fig2.set_size_inches(8, 8)

Writer = ani.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

anim = ani.FuncAnimation(fig2, update2, 
                               frames = int(Nt_r/A), 
                               blit = False, interval=1/30)

anim.save('phase_space_cont.mp4', writer=writer, dpi=100)
    

    
    
   