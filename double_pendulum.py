# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 11:55:40 2021

@author: llucv
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

m=1 #mass of the two bodies
L=1 #lenght of the two 0 mass strings 
g=10 #gravity

def R_K_4(y0,t0,h,f):
    k1=f(t0,y0)
    k2=f(t0+h/2,y0+k1*h/2)
    k3=f(t0+h/2,y0+k2*h/2)
    k4=f(t0+h,y0+k3*h)
    y1=y0+h*(k1+2*k2+2*k3+k4)/6
    return y1

def p_1_d(w1,w2,a1,a2):
    der_p1=m*L*(2*g*np.cos(a1)-L*w1*w2*np.sin(a1-a2))
    return der_p1

def p_2_d(w1,w2,a1,a2):
    der_p2=m*L*(g*np.cos(a2)+L*w1*w2*np.sin(a1-a2))
    return der_p2

def w_1_L(p1,p2,a1,a2):
    w1=(p1-np.cos(a1-a2)*p2)/(m*L*L*(2-np.cos(a1-a2)*np.cos(a1-a2)))
    return w1

def w_2_L(p1,p2,a1,a2):
    w2=(2*p2-np.cos(a1-a2)*p1)/(m*L*L*(2-np.cos(a1-a2)*np.cos(a1-a2)))
    return w2

def d_w_1(a1,a2,w1,w2):
    A=a1-a2
    B=a1-2*a2
    num=-g*(3*np.sin(a1)+np.sin(B))-2*L*np.sin(A)*(w2*w2+w1*w1*np.cos(A))
    den=L*(3-np.cos(2*A))
    d_w1=num/den
    return d_w1

def d_w_2(a1,a2,w1,w2):
    A=a1-a2
    num=2*np.sin(A)*(2*w1*w1*L+2*g*np.cos(a1)+w2*w2*L*np.cos(A))
    den=L*(3-np.cos(2*A))
    d_w2=num/den
    return d_w2

def w_1(a1,a2,w1,w2):
    w_1=w1
    return w_1

def w_2(a1,a2,w1,w2):
    w_2=w2
    return w_2

def der_f(y):
    f=np.zeros((4))
    a1=y[0]
    a2=y[1]
    w1=y[2]
    w2=y[3]
    f[0]=w1
    f[1]=w2
    f[2]=d_w_1(a1,a2,w1,w2)
    f[3]=d_w_2(a1,a2,w1,w2)
    return f 

def R_K_dbl_pnd(y,h):
    #def y: y[0]=a1,y[1]=a2,y[2]=w1,y[3]=w2
    k1=der_f(y)
    k2=der_f(y+h*k1/2)
    k3=der_f(y+h*k2/2)
    k4=der_f(y+k3*h)
    
    y1=y+h*(k1+2*k2+2*k3+k4)/6
    return y1

def double_pendulum(Nt,h,a1_0,a2_0,w1_0,w2_0):
    dbl_pnd=np.zeros((Nt,5))
    dbl_pnd[:,0]=np.arange(Nt)*h
    dbl_pnd[0,1]=a1_0
    dbl_pnd[0,2]=a2_0
    dbl_pnd[0,3]=w1_0
    dbl_pnd[0,4]=w2_0

    
    y=np.zeros((4))
    y[:]=dbl_pnd[0,1:]
      
    for i in range(Nt-1):
        y=R_K_dbl_pnd(y, h)
        a1=y[0]
        a2=y[1]
        w1=y[2]
        w2=y[3]
        Ene=energy(a1, a2, w1, w2)
        print(Ene[0])
        dbl_pnd[i+1,1:]=y[:]
        
    return dbl_pnd
        
def energy(a1,a2,w1,w2):
    Ene=np.zeros((3))
    Ene[1]=m*L*L*(w1*w1+w2*w2*0.5+w1*w2*np.cos(a1-a2))
    Ene[2]=-m*g*L*(2*np.cos(a1)+np.cos(a2))
    Ene[0]=Ene[1]+Ene[2]
    
    return Ene

Nt=2000
h=0.01
a1_0=np.pi/1.5
a2_0=np.pi/4-0.3
w1_0=0
w2_0=0

dbl_pnd=double_pendulum(Nt, h, a1_0, a2_0, w1_0, w2_0)

dbl_pnd_xy=np.zeros((Nt,4))
dbl_pnd_xy[:,0]=L*np.sin(dbl_pnd[:,1])
dbl_pnd_xy[:,1]=-L*np.cos(dbl_pnd[:,1])
dbl_pnd_xy[:,2]=L*np.sin(dbl_pnd[:,2])+L*np.sin(dbl_pnd[:,1])
dbl_pnd_xy[:,3]=-L*np.cos(dbl_pnd[:,2])-L*np.cos(dbl_pnd[:,1])

dbl_pnd_xyt=np.zeros((4))
dbl_pnd_traj=np.zeros((Nt,4))

def update(frame):
    k=frame*5
    dbl_pnd_xyt[:]=dbl_pnd_xy[k,:]
    plt.cla()
    plt.plot(dbl_pnd_xy[0:k,0],dbl_pnd_xy[0:k,1],'r-')
    plt.plot(dbl_pnd_xy[0:k,2],dbl_pnd_xy[0:k,3],'b-')
    
    x_line_12 = [dbl_pnd_xyt[0], dbl_pnd_xyt[2]]
    y_line_12 = [dbl_pnd_xyt[1], dbl_pnd_xyt[3]]

    x_line_01 = [dbl_pnd_xyt[0], 0]
    y_line_02 = [dbl_pnd_xyt[1], 0]

    plt.plot(x_line_12, y_line_12)
    plt.plot(x_line_01, y_line_02)
    plt.plot(dbl_pnd_xyt[0],dbl_pnd_xyt[1],'ro')
    plt.plot(dbl_pnd_xyt[2],dbl_pnd_xyt[3],'bo')
    plt.axis([-2.2,2.2,-2.2,2.2])

fig = plt.figure(figsize=(8, 8))
ax1 = plt.subplot()

Writer = ani.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

anim = ani.FuncAnimation(fig, update, 
                               frames = int((Nt-1)/5), 
                               blit = False, interval=100)

anim.save(str(a1_0)+"_"+str(a2_0)+'double_pendulum.mp4', writer=writer)



    

