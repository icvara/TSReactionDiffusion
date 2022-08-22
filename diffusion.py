'''
Here are the script to run stuff for diffusion
written by Içvara


diffusion resolve with explicit method of Forward Euler scheme 
the issue with this methods is to not have error and good precision the dt need to be very small or the dx very big
I try to solve it by using anoth layer with lower dx precision


https://hplgit.github.io/fdm-book/doc/pub/book/sphinx/._book011.html
'''
import numpy as np
import matplotlib.pyplot as plt

dx_diffusion=0.5

def laplacian1D(a, dx):
    return (
        - 2 * a
        + np.roll(a,1,axis=0) 
        + np.roll(a,-1,axis=0)
    ) / (dx ** 2)


def laplacian2D(a, dx):
    return (
        - 4 * a
        + np.roll(a,1,axis=0) 
        + np.roll(a,-1,axis=0)
        + np.roll(a,+1,axis=1)
        + np.roll(a,-1,axis=1)
    ) / (dx ** 2)



def diffusion1D(ui,D,dx,dt):
	ud = laplacian1D(ui, dx)
	return dt * (D * ud)

def diffusion2D(ui,D,dx,dt):
    ud = laplacian2D(ui, dx)
    return dt * (D * ud)

def increase_matrix(A0,dx,dx_diffusion):
    time=int(dx_diffusion/dx)
    A0_extended_int=np.repeat(A0,time,axis=1)
    A0_extended=np.repeat(A0_extended_int,time,axis=0)
    return A0_extended

def decrease_matrix(A0,dx,dx_diffusion):
    time=int(dx_diffusion/dx)
    A0_decreased = A0[::time,::time]
    return A0_decreased

'''
L=5
u0=np.ones((L,L))
u0[2,2]=10

plt.imshow(u0)
plt.show()
tt=50
utot=np.zeros((tt,L,L))
utot[0]=u0
u=u0
for i in np.arange(tt):
    u = u  + diffusion2D(u,D=2)
    utot[i] = u 



plt.plot(utot[:,2,2])
plt.show()

print(utot)
for i in np.arange(tt):
    plt.plot(utot[i,:,2])
plt.show()
'''