import numpy as np
import matplotlib.pyplot as plt
import random

'''
width=50
dt=0.001
dx=0.1

time=0.1

'''

par_growth={  
    'D_growth': 1, #density growing by hour with consideration of the max density
    'max_density': 100,
	'H_growth': 0.01 #  cm/h ?
}

def calculateStationary(d0,par):
	stat=np.copy(d0)
	'''
	stat[stat<=0.9*par['max_density']]=0
	stat[stat>0.9*par['max_density']]=1
	'''
	stat = stat/par['max_density']
	return stat

def densityGrowth_model(d0,par):
	d =  par['D_growth']*d0 * (1 - d0 / par['max_density'])
	return d

def horizontalGrowth_model(c0,center,t,dx,par):
	r= par['H_growth']*t /dx
	cell=c0.copy()
	init=False
	for i in np.arange(len(center[0])):	
		centerx=center[0][i]
		centery=center[1][i]
		emptycell = np.argwhere(c0==0)
		dist = np.sqrt(np.power(emptycell[:,0]-centerx,2)+np.power(emptycell[:,1]-centery,2))
		insidecircle =  np.argwhere(dist<r)
		if init==False:
			gcell=emptycell[insidecircle]
			init=True
		else:
			gcell = np.vstack([gcell, emptycell[insidecircle]])
	return gcell

def roll(C0,dist):
	left=np.roll(C0,dist,axis=1)
	right=np.roll(C0,-dist,axis=1)
	down=np.roll(C0,-dist,axis=0)
	up=np.roll(C0,dist,axis=0)
	leftbottom=np.roll(left,-dist,axis=0)
	rightbottom=np.roll(right,-dist,axis=0)
	lefttop=np.roll(left,dist,axis=0)
	righttop=np.roll(right,dist,axis=0)

	return left,right, down, up, leftbottom, rightbottom, lefttop, righttop


def heredity(G0,R0,A0,C0,gcell,size,ngh=1):

	left,right, down, up, leftbottom, rightbottom, lefttop, righttop =  roll(C0,ngh)
	Gleft,Gright, Gdown, Gup, Gleftbottom, Grightbottom, Glefttop, Grighttop =  roll(G0,ngh)
	Rleft,Rright, Rdown, Rup, Rleftbottom, Rrightbottom, Rlefttop, Rrighttop =  roll(R0,ngh)
	Aleft,Aright, Adown, Aup, Aleftbottom, Arightbottom, Alefttop, Arighttop =  roll(A0,ngh)




	hg = np.zeros((size,size))
	hr = np.zeros((size,size))
	ha = np.zeros((size,size))



	for gri,grc in enumerate(gcell[:,0,1]):
			i = gcell[gri,:,0]
			j = gcell[gri,:,1]
		
			#index_array=np.array(([i+1,j],[i-1,j],[i,j+1],[i,j-1],[i+1,j+1],[i-1,j-1],[i+1,j-1],[i-1,j+1]))
			C0_array=np.array((left[i,j],right[i,j],down[i,j],up[i,j],leftbottom[i,j],rightbottom[i,j],lefttop[i,j],righttop[i,j]))
			G0_array=np.array((Gleft[i,j],Gright[i,j],Gdown[i,j],Gup[i,j],Gleftbottom[i,j],Grightbottom[i,j],Glefttop[i,j],Grighttop[i,j]))
			R0_array=np.array((Rleft[i,j],Rright[i,j],Rdown[i,j],Rup[i,j],Rleftbottom[i,j],Rrightbottom[i,j],Rlefttop[i,j],Rrighttop[i,j]))
			A0_array=np.array((Aleft[i,j],Aright[i,j],Adown[i,j],Aup[i,j],Aleftbottom[i,j],Arightbottom[i,j],Alefttop[i,j],Arighttop[i,j]))

			thereisacell = np.where(C0_array==1)
			if len(thereisacell[0])>0:
				index = random.choice(thereisacell[0])
				hg[i,j]=G0_array[index]
				hr[i,j]=R0_array[index]
				ha[i,j]=A0_array[index]


			else:
				ngh += 1
				left,right, down, up, leftbottom, rightbottom, lefttop, righttop =  roll(C0,ngh)
				Gleft,Gright, Gdown, Gup, Gleftbottom, Grightbottom, Glefttop, Grighttop =  roll(G0,ngh)
				Rleft,Rright, Rdown, Rup, Rleftbottom, Rrightbottom, Rlefttop, Rrighttop =  roll(R0,ngh)
				Aleft,Aright, Adown, Aup, Aleftbottom, Arightbottom, Alefttop, Arighttop =  roll(A0,ngh)


				C0_array=np.array((left[i,j],right[i,j],down[i,j],up[i,j],leftbottom[i,j],rightbottom[i,j],lefttop[i,j],righttop[i,j]))
				G0_array=np.array((Gleft[i,j],Gright[i,j],Gdown[i,j],Gup[i,j],Gleftbottom[i,j],Grightbottom[i,j],Glefttop[i,j],Grighttop[i,j]))
				R0_array=np.array((Rleft[i,j],Rright[i,j],Rdown[i,j],Rup[i,j],Rleftbottom[i,j],Rrightbottom[i,j],Rlefttop[i,j],Rrighttop[i,j]))
				A0_array=np.array((Aleft[i,j],Aright[i,j],Adown[i,j],Aup[i,j],Aleftbottom[i,j],Arightbottom[i,j],Alefttop[i,j],Arighttop[i,j]))


				thereisacell = np.where(C0_array==1)
				if len(thereisacell[0])>0:
						index = random.choice(thereisacell[0])
						hg[i,j]=G0_array[index]
						hr[i,j]=R0_array[index]
						ha[i,j]=A0_array[index]

				else:
						print("yolo: growth too fast no neighbourg cell")



	return hg,hr,ha







'''
nx=round(width/dx)
d0=np.zeros((nx,nx))
c0=np.zeros((nx,nx))
G0=np.zeros((nx,nx))
R0=np.zeros((nx,nx))
#d0[200,249]=1
#c0[20,249]=1

c0[round(nx/2),round(nx/2)]=1
G0[round(nx/2),round(nx/2)]=1

c0[round(nx/2)+2,round(nx/2)-2]=1
R0[round(nx/2)+2,round(nx/2)-2]=1

c0[round(nx/2)+2,round(nx/2)+2]=1


c0[nx-1,round(nx/2)+2]=1






par=par_growth

C,D,G,R =simulation(c0,d0,G0,R0,par,width,dx,time,dt)

fig,axs= plt.subplots(3,4,sharey='none')
axs[0,0].imshow(C[0])
axs[1,0].imshow(G[0],cmap="Greens")
axs[2,0].imshow(R[0],cmap="Reds")

axs[0,1].imshow(C[5])
axs[1,1].imshow(G[5],cmap="Greens")
axs[2,1].imshow(R[5],cmap="Reds")

axs[0,2].imshow(C[10])
axs[1,2].imshow(G[10],cmap="Greens")
axs[2,2].imshow(R[10],cmap="Reds")

axs[0,3].imshow(C[-1])
axs[1,3].imshow(G[-1],cmap="Greens")
axs[2,3].imshow(R[-1],cmap="Reds")

plt.show()

plt.show()

plt.plot(D[-1,249,:])
plt.show()



'''

