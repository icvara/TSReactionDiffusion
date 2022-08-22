import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import seaborn as sns


par={

	'beta_green':1,
	'K_ahl_green':1.5,
	'n_ahl_green':2,
	'delta_green':1,

	#'K_GREEN':1000,
	#'n_GREEN':2,

	'beta_red':1.00,
	'K_ahl_red':1,
	'n_ahl_red':2,
	'delta_red':1.,

	'K_RED':2.5,
	'n_RED':2,


    'D_red':1,
    'D_green':0.01

    }

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




def model_turing(GREENi,REDi,par):

    GREEN = 0.01 + par['beta_green']*np.power(GREENi*10**par['K_ahl_green'],par['n_ahl_green'])/(1+np.power(GREENi*10**par['K_ahl_green'],par['n_ahl_green']))
    GREEN = GREEN / (1 + np.power(REDi*10**par['K_RED'],par['n_RED']))
    GREEN = GREEN - par['delta_green']*GREENi  
    

    RED = (par['beta_red']*np.power(GREENi*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(GREENi*10**par['K_ahl_red'],par['n_ahl_red']))
   # RED = RED / (1 + np.power(GREENi*10**par['K_GREEN'],par['n_GREEN']))
    RED = RED - par['delta_red']*REDi 

    return GREEN,RED




def simulation(Gi,Ri,par,dt,tt=1):
	ti=1
	av=np.zeros((round(tt/dt),2,len(Gi)))
	av[0,0]=Gi
	av[0,1]=Ri


	while (ti < round(tt/dt)):
			G,R =  model_turing(Gi,Ri,par)
			Rd = laplacian1D(Ri, dx)
			Gd = laplacian1D(Gi, dx)

			#delta_a = dt * (D * La)       
			#a += delta_a
			#at[ti]=a

			Gi += dt*G + dt * (par['D_green'] * Gd)
			Ri += dt*R + dt * (par['D_red'] * Rd)
			av[ti,0]=Gi 
			av[ti,1]=Ri 

			ti += 1

	return av


def simulation2D(Gi,Ri,par,dt,tt=1):
	ti=1
	#av=np.zeros((int(tt/dt),2,len(Gi[0]),len(Gi[1])))
	#av[0,0]=Gi
	#av[0,1]=Ri
	

	while (ti < int(tt/dt)):
			G,R =  model_turing(Gi,Ri,par)
			Rd = laplacian2D(Ri, dx)
			Gd = laplacian2D(Gi, dx)

			#delta_a = dt * (D * La)       
			#a += delta_a
			#at[ti]=a

			Gi += dt*G + dt * (par['D_green'] * Gd)
			Ri += dt*R + dt * (par['D_red'] * Rd)

			#av[ti,0]=Gi 
			#av[ti,1]=Ri 

			ti += 1
	av=np.array((Gi,Ri))
	print(av.shape)

	return av




width=200
dx=10/width
D=1
dt= 0.9 * (dx ** 2) / (2 * D)
print(dt)

F = D*dt/dx**2  
print("F values is: " + str(F))

time=100
G=np.zeros(width)
G[0:40]=0.2
G=np.random.normal(loc=0, scale=0.5, size=width)
R=np.zeros(width)
R[50:-1]=0.2
R=np.random.normal(loc=0, scale=0.5, size=width)
av=simulation(G,R,par,dt,tt=time)

fig, axs = plt.subplots(5)
tt=int(time/dt/4)
print(tt)

axs[0].plot(av[0,0],'--g')
axs[0].plot(av[0,1],'--r')
axs[1].plot(av[tt,0],'g')
axs[1].plot(av[tt,1],'r')
axs[2].plot(av[(tt*2),0],'g')
axs[2].plot(av[(tt*2),1],'r')
axs[3].plot(av[(tt*3),0],'g')
axs[3].plot(av[(tt*3),1],'r')
axs[4].plot(av[-1,0],'g')
axs[4].plot(av[-1,1],'r')
#plt.plot(at[round(tt/2)])
#plt.plot(at[-1])
plt.show()


G=np.ones((width,width))
R=np.ones((width,width))
for i in np.arange(len(R[0])):
	G[i,:]=G[i,:]*np.random.normal(loc=0, scale=0.5, size=width)
	R[i,:]=R[i,:]*np.random.normal(loc=0, scale=0.5, size=width)



av=simulation2D(G,R,par,dt,tt=time)

fig, axs = plt.subplots(2,2)

tt=int(time/dt/4)

sns.heatmap(ax=axs[0,0], data= av[0])


sns.heatmap(ax=axs[0,1], data= av[1])


#plt.imshow(av[0,0])
'''
axs[0].plot(av[0,1],'--r')
axs[1].plot(av[tt,0],'g')
axs[1].plot(av[tt,1],'r')
axs[2].plot(av[(tt*2),0],'g')
axs[2].plot(av[(tt*2),1],'r')
axs[3].plot(av[(tt*3),0],'g')
axs[3].plot(av[(tt*3),1],'r')
axs[4].plot(av[-1,0],'g')
axs[4].plot(av[-1,1],'r')
'''
#plt.plot(at[round(tt/2)])
#plt.plot(at[-1])
plt.show()



'''
def animate(i):
	y=av[i,0]
	line.set_ydata(y)
	return line,

fig, ax = plt.subplots()
line, = ax.plot(av[0,0])

ani = animation.FuncAnimation(
    fig, animate, interval=0.00001, blit=True)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# writer = animation.FFMpegWriter(
#     fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

plt.show()
'''

