'''
================================================================================================================================================================

Script to make diverse plot

==============================================================================================================================================================
'''

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LogNorm

import numpy as np
import os

from model import *
from grow import *
from torun import *

import imageio
import seaborn as sns
import pandas as pd



#coltor blind color
color=sns.color_palette("colorblind")
colorGREEN=color[2]
colorBLUE=color[0]
colorRED=color[3]
colorORANGE=color[1]
colorCYAN=color[9]


colorPurple=color[4]
'''
for i in np.arange(len(color)):
	print(color[i])
	plt.plot(i,i,'o',color=color[i])
plt.show()
'''






def generate_gif(av):
#this doesn't work so well
	fig, axs = plt.subplots(2,5)
	filenames = []
	frames = []

	for i in np.arange(av.shape[0]):
		axs = plot2D_simple(av,-1, axs, 1, 0)
		axs = plot2D_simple(av,i, axs, 0, 0)
		filename = f'{i}.png'
		filenames.append(filename)
		fig.savefig(filename)

	for filename in filenames:
	        frames.append(imageio.imread(filename))


	# Save them as frames into a gif 
	exportname = "output.gif"
	kargs = { 'duration': 0.001 }
	imageio.mimsave(exportname, frames, 'GIF', **kargs)       
# Remove files
	for filename in set(filenames):
	    os.remove(filename)


'''
=================================================================================================================================================================================

2D grid plot

=================================================================================================================================================================================
'''

def bestMatrix(G,R):

	Mm=np.ones(G.shape)*0#np.nan
	#Rm=np.ones(R.shape)*0#np.nan

	for i in np.arange(G.shape[0]):
		for j in np.arange(G.shape[1]):
			if np.isnan(G[i,j])==False:
				if G[i,j]>R[i,j]:
					Mm[i,j]=G[i,j]
				elif G[i,j]<R[i,j]:
					Mm[i,j]=R[i,j]*-1
				elif G[i,j]==R[i,j]:
					Mm[i,j]=0#R[i,j]



	return Mm

def bestMatrix3(G,R,B):

	Gm=np.ones(G.shape)*np.nan
	Rm=np.ones(R.shape)*np.nan
	Bm=np.ones(B.shape)*np.nan



	idxG=np.where(G>B)
	idxR=np.where(R>G)
	idxB=np.where(B>R)



	Gm[idxG]=G[idxG]
	Rm[idxR]=R[idxR]
	Bm[idxB]=B[idxB]

	return Gm,Rm,Bm

def plot2D_simple(av,tt, axs, i, j, n=5):
	colorlist=['Greens','Reds','Blues','Greys','viridis']
	maxV=1
	for ni in np.arange(n):
			if ni == 4:
				maxV=100
			axs[i,j+ni].imshow(av[tt,ni],cmap=colorlist[ni],vmin=0,vmax=maxV)

	return axs

def plot2D_kymograph(av,w, axs, i, j,n=1):
	maxV=1
	colorlist=['Greens','Reds','Blues','Greys','viridis']
	for ni in np.arange(n):
			if ni == 4:
				maxV=100
			axs[i,j+ni].imshow(av[:,ni,:,w].T,cmap=colorlist[ni],aspect='auto')#,vmin=0,vmax=maxV)
	return axs


def plot2D_simple_tgh(av,tt, axs, i, j):

	colorlist=['Greens','Reds','Blues','Greys','viridis']
	Mm = bestMatrix (av[tt,0],av[tt,1])
	NN=Mm.copy()
	NN[NN!=0]=np.nan


	axs[i,j].imshow(Mm,cmap="RdYlGn",alpha=1,vmin=-100,vmax=100)#,vmin=0,vmax=vmaxG)

	#axs[i,j].imshow(NN,cmap="Greys",alpha=1)#,vmin=0,vmax=vmaxG)

	#axs[i,j].imshow(Rm,cmap=colorlist[1],alpha=0.5)#,vmin=0,vmax=vmaxG)


	return axs


def plot2D_simple_tgh3(av,tt, axs, i, j):

	colorlist=['Greens','Reds','Blues','Greys','viridis']
	Gm,Rm,Bm = bestMatrix3 (av[tt,0],av[tt,1],av[tt,3])
	vmaxG=np.nanmax(Gm)
	vmaxR=np.nanmax(Rm)
	vmaxB=np.nanmax(Bm)
	
	axs[i,j].imshow(Bm,cmap=colorlist[2],alpha=1,vmin=0,vmax=1)#vmin=0,vmax=1
	axs[i,j].imshow(Gm,cmap=colorlist[0],alpha=1,vmin=0,vmax=1)
	axs[i,j].imshow(Rm,cmap=colorlist[1],alpha=1,vmin=0,vmax=1)#vmin=0,vmax=1
	
	#axs[i,j].imshow(av[tt,3],cmap=colorlist[2],alpha=1)#vmin=0,vmax=vmaxB)#vmin=0,vmax=1
	#axs[i,j].imshow(av[tt,0],cmap=colorlist[0],alpha=0.5)#,vmin=0,vmax=vmaxG)
	#axs[i,j].imshow(av[tt,1],cmap=colorlist[1],alpha=0.5)#vmin=0,vmax=vmaxR)#vmin=0,vmax=1

	return axs

def plot2D_kymograph_tgh(av,w, axs, i, j):
	Mm = bestMatrix (av[:,0,:,w],av[:,1,:,w])



	colorlist=['Greens','Reds','Blues','Greys','viridis']
	axs[i,j].imshow(Mm.T,cmap="RdYlGn",alpha=1,aspect='auto')
	#axs[i,j].imshow(Rm.T,cmap=colorlist[1],alpha=1,aspect='auto')
	return axs


def plot_crossection(av,tt,w, axs, i, j):

	axs[i,j].plot(av[tt,4,:,w],'k')

	axs[i,j].plot(av[tt,0,:,w],color=colorGREEN)
	axs[i,j].plot(av[tt,1,:,w],color=colorRED)
	#axs[i,j].plot(av[tt,2,:,w],color=colorORANGE)

	axs[i,j].plot(av[tt,3,:,w],color=colorBLUE)

	axs[i,j].plot(av[tt,5,:,w],'m')

	axs[i,j].set_yscale('log')

	#axs[i,j].plot(av[tt,4,:,w])
	return axs

def plot_crossection_diffusion(av,tt,w, axs, i, j):
	axs[i,j].plot(av[tt,4,:,w],'k')
	axs[i,j].plot(av[tt,3,:,w],color=colorBLUE)
	#axs[i,j].plot(av[tt,5,:,w]/np.max(av[tt,5,:,w]),'m')
	#axs[i,j].plot(av[tt,1,:,w]/np.max(av[tt,1,:,w]),color=colorGREEN)
	#axs[i,j].plot(av[tt,2,:,w]/np.max(av[tt,2,:,w]),color=colorRED)
	#axs[i,j].plot(av[tt,3,:,w]/np.max(av[tt,3,:,w]),color=colorORANGE)

	axs[i,j].plot(av[tt,5,:,w],'m')
	axs[i,j].plot(av[tt,1,:,w],color=colorGREEN)
	axs[i,j].plot(av[tt,2,:,w],color=colorRED)
	axs[i,j].plot(av[tt,3,:,w],color=colorORANGE)


	#axs[i,j].plot(av[tt,4,:,w])
	return axs

def plot_crosstime(av,w, axs, i, j):

	axs[i,j].plot(av[:,0,w,w],color=colorGREEN)
	axs[i,j].plot(av[:,1,w,w],color=colorRED)
	axs[i,j].plot(av[:,2,w,w],color=colorORANGE)

	axs[i,j].plot(av[:,3,w,w],color=colorBLUE)
	#axs[i,j].plot(av[:,5,w,w]/np.max(av[:,5,w,w]),color='m')

	#axs[i,j].plot(av[:,4,w,w],'k')
	#axs[i,j].plot(av[:,4,w,w])
	return axs


def test_diffusion(name):
	p={
		'beta_green':0,
		'alpha_green':1,#0.001,

		'K_ahl_green':0,
		'n_ahl_green':2,
		'delta_green':1,

		'K_GREEN':1,
		'n_GREEN':2,

		'beta_red':0,
		'alpha_red':0.001,

		'K_ahl_red':0,
		'n_ahl_red':2,
		'delta_red':1.,

		'K_RED':-10,
		'n_RED':2,

		'beta_ahl':1,
		'K_ahl':1,
		'n_ahl':2,
		'delta_ahl':1,

	    'D_ahl':1,#0.1,

	    'K_IPTG':1
    }

	time=100#20
	dt=0.005
	size=4
	dx=0.05
	width=int(size/dx)
	pos=int(width/2)


	par_growth['H_growth']=0

	G=np.ones((width,width))*0
	R=np.ones((width,width))*0
	A=np.ones((width,width))*0
	C=np.ones((width,width))*0
	D=np.ones((width,width))*0
	A0=np.zeros(1)
	IPTG=np.zeros(1)
	fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(12,4))

	#colist=['k','r','g','b','m']


	k=np.logspace(-1,-2,4)
	k=[0.1]
	#k= 5*10**-6 *3600  #cm2/s -> cm2/h	4.9 x 10‑6 cm2/s	AHL diffusion constant	Stewart P.S., 2003 , this values seems to slow
	#k=np.linspace(0,1,5)
	nx=int(1/dx)
	x1= np.arange(1,nx)*2 *dx
	for c,Dif in enumerate(k):
		center=[]
		border=[]
		AHLdist=[]
		p['D_ahl']=Dif
		for x in np.arange(1,nx):
			#for i in np.arange(-x,x+1):
			#	for j in np.arange(-x,x+1):
			C=np.ones((width,width))*0
			G=np.ones((width,width))*0
			R=np.ones((width,width))*0
			A=np.ones((width,width))*0
			C=np.ones((width,width))*0
			D=np.ones((width,width))*0


			index=drawcircle(x,C,[pos,pos])
			C[index[:,0,0],index[:,0,1]]=1
			D[index[:,0,0],index[:,0,1]]=1
			G[index[:,0,0],index[:,0,1]]=1
			R[index[:,0,0],index[:,0,1]]=0.01

					#C[pos+i,pos+j]=1
					#D[pos+i,pos+j]=1
					#G[pos+i,pos+j]=1
					#R[pos+i,pos+j]=0.01

			av=simulation(G,R,A,C,D,A0,IPTG,p,width,dx,time,dt,"TSXLT",time,False,False)
			center.append(av[-1,3,pos,pos])
			border.append(av[-1,3,pos,pos]-av[-1,3,pos+1-x,pos])
			#AHLdist.append(av[-1,3,pos,:])
		#for l in AHLdist:
			axs[2].plot(av[-1,4,pos,:],color=color[c])

			axs[2].plot(av[-1,3,pos,:],color=color[c])



		axs[0].plot(x1,center,'-o',color=color[c])
		axs[1].plot(x1,border,'-o',color=color[c])



	axs[0].set_xlabel("Colony size [cm]")
	axs[1].set_title("duration[h]" + str(time) + " with " + str(dt)  + " dt " + str(dx) + " dx")
	axs[0].set_ylabel("AHL at center")
	axs[1].set_xlabel("Colony size [cm]")
	axs[1].set_ylabel("AHL center - AHL border")
	axs[1].legend( k, loc=0)
	plt.savefig("Figures/"+name+"_testDiffusion.pdf",dpi=300)
	plt.show()


def drawcircle(r,c0,center):
	centerx=center[0]
	centery=center[1]
	emptycell = np.argwhere(c0==0)
	dist = np.sqrt(np.power(emptycell[:,0]-centerx,2)+np.power(emptycell[:,1]-centery,2))
	insidecircle =  np.argwhere(dist<r)
	idx=emptycell[insidecircle]
	#c0[insidecircle]=1
	return idx


'''
=================================================================================================================================================================================

bifurcation plot

=================================================================================================================================================================================
'''



def bifu_plot(par,A0,axs,i,j):
	#A0[0]=0
	colors = [colorBLUE, colorGREEN] # first color is black, last is red


	ss=findss(A0,[0],par)
	sse=getEigen(ss,A0,par)
	I=0

	for ai,a in enumerate(A0):
		for si in np.arange(0,ss.shape[2]):
			e=sse[ai,I,si,:]
			er=e.real
			ei=e.imag
			

			sss=ss[ai,I,si,:]
			if np.all(er!=0): 
				#print(sse)
				if np.all(er<0): #stable
					#if np.all(ei==0):
					#node
						m='o'
						axs[i,j].plot(np.log10(a),sss[0],m, c=colorGREEN ,markersize=1)
						axs[i,j].plot(np.log10(a),sss[1],m, c=colorRED, markersize=1.)
						#axs[i,j+2].plot(np.log10(a),sss[2],m, c= colorBLUE ,markersize=1.)


				elif np.any(er>0): #unstable
					#check complex conjugate
					if len(er) != len(set(er)):
						if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
							#oscilation
							m='o'
							axs[i,j+0].plot(np.log10(a),sss[0],m+'m',markersize=1.)
							axs[i,j].plot(np.log10(a),sss[1],m+'m',markersize=1.)
							#axs[i,j+2].plot(np.log10(a),sss[2],m+'m',markersize=1.)
						else:
							m='x'
							axs[i,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
							axs[i,j].plot(np.log10(a),sss[1],m+'k',markersize=1.)
							#axs[i,j+2].plot(np.log10(a),sss[2],m+'k',markersize=1.)
					else:
							m='x'
							axs[i,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
							axs[i,j].plot(np.log10(a),sss[1],m+'k',markersize=1.)
							#axs[i,j+2].plot(np.log10(a),sss[2],m+'k',markersize=1.)

	#axs[i,j].set_ylim(-0.1,1.1)
	#axs[i,j+1].set_ylim(-0.1,1.1)
	#axs[i,j+2].set_ylim(-0.1,1.1)

	return axs

def bifu_plot_tgh(par,A0,axs,i,j):
	ss=findss(A0,[0],par)
	sse=getEigen(ss,A0,par)
	I=0
	for ai,a in enumerate(A0):
		for si in np.arange(0,ss.shape[2]):
			e=sse[ai,I,si,:]
			er=e.real
			ei=e.imag
			
			sss=ss[ai,I,si,:]
			if np.all(er!=0): 
				#print(sse)
				if np.all(er<0): #stable
					#if np.all(ei==0):
					#node
						m='o'
						axs[i,j].plot(np.log10(a),sss[0],m, c=colorGREEN ,markersize=2)
						axs[i,j].plot(np.log10(a),sss[1],m, c=colorRED, markersize=2.)
						axs[i,j].plot(np.log10(a),sss[2],m, c= colorBLUE ,markersize=2.)
				elif np.any(er>0): #unstable
					#check complex conjugate
					if len(er) != len(set(er)):
						if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
							#oscilation
							m='o'
							axs[i,j].plot(np.log10(a),sss[0],m,c=colorGREEN,markersize=2., mfc='none')
							axs[i,j].plot(np.log10(a),sss[1],m,c=colorRED,markersize=2., mfc='none')
							axs[i,j].plot(np.log10(a),sss[2],m,c=colorBLUE,markersize=2., mfc='none')
							axs[i,j].plot(np.log10(a),sss[0],m,c='m',markersize=.5)
							axs[i,j].plot(np.log10(a),sss[1],m,c='m',markersize=.5)
							axs[i,j].plot(np.log10(a),sss[2],m,c='m',markersize=.5)
						else:
							m='x'
							axs[i,j].plot(np.log10(a),sss[0],m,c=colorGREEN,markersize=.5)
							axs[i,j].plot(np.log10(a),sss[1],m, c=colorRED,markersize=.5)
							axs[i,j].plot(np.log10(a),sss[2],m, c=colorBLUE,markersize=.5)
					else:
							m='x'
							axs[i,j].plot(np.log10(a),sss[0],m,c=colorGREEN,markersize=.5)
							axs[i,j].plot(np.log10(a),sss[1],m, c=colorRED,markersize=.5)
							axs[i,j].plot(np.log10(a),sss[2],m, c=colorBLUE,markersize=.5)

	return axs
def bifu_plot_par(par,A0,axs,i,j,k,s):
	#A0[0]=0
	custom_colors = [colorRED, colorGREEN] # first color is black, last is red
	col = colors.LinearSegmentedColormap.from_list("Custom", custom_colors)
	col=cm.viridis
	color=np.linspace(0, 1,num=s)
	ci=color[k]

	ss=findss(A0,[0],par)
	sse=getEigen(ss,A0,par)
	I=0

	for ai,a in enumerate(A0):
		for si in np.arange(0,ss.shape[2]):
			e=sse[ai,I,si,:]
			er=e.real
			ei=e.imag
			

			sss=ss[ai,I,si,:]
			if np.all(er!=0): 
				#print(sse)
				if np.all(er<0): #stable
					#if np.all(ei==0):
					#node
						m='o'
						axs[i,j].plot(np.log10(a),sss[0],m,c=col(ci),markersize=1)
						axs[i,j].plot(np.log10(a),sss[1],m,c='r',markersize=1)

						#axs[i+1,j].plot(np.log10(a),sss[1],m+'r',markersize=1.)
						#axs[i+2,j].plot(np.log10(a),sss[2],m+'b',markersize=1.)


				elif np.any(er>0): #unstable
					#check complex conjugate
					if len(er) != len(set(er)):
						if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
							#oscilation
							m='o'
							axs[i,j+0].plot(np.log10(a),sss[0],m, c=col(ci),markersize=4., mfc='none')
							#axs[i+1,j].plot(np.log10(a),sss[1],m+'m',markersize=1.)
							#axs[i+2,j].plot(np.log10(a),sss[2],m+'m',markersize=1.)
						else:
							m='x'
							axs[i,j+0].plot(np.log10(a),sss[0],m, c=col(ci),markersize=.5)
							#axs[i+1,j].plot(np.log10(a),sss[1],m+'k',markersize=1.)
							#axs[i+2,j].plot(np.log10(a),sss[2],m+'k',markersize=1.)
					else:
							m='x'
							axs[i,j+0].plot(np.log10(a),sss[0],m, c=col(ci),markersize=.5)
							#axs[i+1,j].plot(np.log10(a),sss[1],m+'k',markersize=1.)
							#axs[i+2,j].plot(np.log10(a),sss[2],m+'k',markersize=1.)

	#axs[i,j].set_ylim(-0.1,1.1)
	#axs[i,j+1].set_ylim(-0.1,1.1)
	#axs[i,j+2].set_ylim(-0.1,1.1)

	return axs



class MidpointNormalize(mpl.colors.Normalize):
	#class to set up the midpoint of colorbar
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
    	x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
    	return np.ma.masked_array(np.interp(value, x, y))


def isOscillation(e):
	b=0
	if np.any(e.real>0): #unstable
					#check complex conjugate
		if len(e.real) != len(set(e.real)):
			if np.any(e.imag==0) and len(e.imag) != len(set(np.abs(e.imag))):
				b=1
	return b


def monostable_plot(A0,par,axs,i,j,k,p):
	mono_matrix=np.zeros((len(A0),len(k)))
	mono_green=np.zeros((len(A0),len(k)))
	mono_red=np.zeros((len(A0),len(k)))
	oscil_matrix=np.zeros((len(A0),len(k)))

	maxvalues=par['beta_green']
	minvalues=0.001#par['alpha_green']


	for ki,kk in enumerate(k):
		par[p]=kk
		ss=findss(A0,[0],par)
		sse=getEigen(ss,A0,par)
		oscil=np.apply_along_axis(isOscillation,3,sse)
		c=np.count_nonzero(~np.isnan(ss[:,:,:,0]),axis=2)
		c_oscil=np.count_nonzero(oscil,axis=2)
		mono_matrix[:,ki]=c[:,0]
		mono_green[:,ki]=ss[:,:,0,0][:,0]#-minvalues)/(maxvalues-minvalues)
		mono_red[:,ki]=ss[:,:,0,1][:,0]#-minvalues)/(maxvalues-minvalues)
		oscil_matrix[:,ki]=c_oscil[:,0]

	mono_matrix[mono_matrix==1]=np.nan
	#mono_green[mono_green<0.05*maxvalues]=np.nan
	#mono_red[mono_red<0.05*maxvalues]=np.nan
	oscil_matrix[oscil_matrix==0]=np.nan

	
	axs[i,j].imshow(mono_red.T,aspect="auto", cmap="Reds",vmin=minvalues,vmax=maxvalues, alpha=1)
	axs[i,j].imshow(mono_green.T,aspect="auto", cmap="Greens",vmin=minvalues,vmax=maxvalues, alpha=0.7)
	axs[i,j].imshow(mono_matrix.T,aspect="auto", cmap="Blues", vmin=0,vmax=5)# vmin=1,vmax=5)
	#axs[i,j].imshow(oscil_matrix.T,aspect="auto", cmap="Greys", vmin=0,vmax=1)# vmin=1,vmax=5)
	a=np.argwhere(oscil_matrix.T==1)
	axs[i,j].plot(a[:,1],a[:,0],'mx',markersize=.5)

	indextick=[]
	for ind in np.arange(5):
		new_i = int(len(k)/5)*ind
		indextick.append(new_i)

	#indextick=np.arange(0,len(k))
	Xlabel=np.linspace(min(k),max(k),5)
	Ylabel=np.logspace(min(A0),max(A0),5)

	axs[i,j].set_yticks(indextick)
	axs[i,j].set_xticks(indextick)
	axs[i,j].set_yticklabels(np.round(k[indextick],2), fontsize=6 )
	axs[i,j].set_xticklabels(np.round(A0[indextick],3),rotation=90,fontsize=6 )

	'''
	axs[i,j].set_yticks(np.arange(len(k)))
	axs[i,j].set_xticks(np.arange(len(A0)))
	axs[i,j].set_yticklabels(np.round(k,2), fontsize=6)
	axs[i,j].set_xticklabels(np.round(A0,2),rotation=90,fontsize=6)
	'''

	return axs

'''
===========================================================================================================================================================================================================================

plot for par space here

=============================================================================================================================================================================================================================
'''



def par_plot(df_par,parlist,cl=None):
    fig, axs = plt.subplots(len(parlist),len(parlist))#, constrained_layout=True, figsize=(10/inch_cm,10/inch_cm))

    for x in np.arange(len(parlist)):
        for y in np.arange(len(parlist)):

            px=parlist[x]['name']
            py=parlist[y]['name']

            if x==y:
                #axs[x,y] = sns.distplot(df_par[px])
                axs[x,y].hist(df_par[px],bins=15, density=True)

                axs[x,y].set_xlim(parlist[x]['lower_limit'],parlist[x]['upper_limit'])
                #sns.kdeplot(ax=axs[x, y],x=df_par[px])

            elif y>x:
                txt=df_par[px][0]
                corr=np.corrcoef(df_par[px],df_par[py])[0, 1] #pearson correlation
                corr_r=np.round(corr,2)
                axs[x,y].imshow([[corr],[corr]], vmin=-1, vmax=1, cmap="bwr",aspect="auto")
                axs[x,y].text(0,0.5,corr_r,fontsize=8,ha='center', va='center')

            elif y<x:
                #z = gaussian_kde(df_par[px])(df_par[py])
                #idx = z.argsort()
                #xx, yy, z = df_par[px][idx], df_par[py][idx], z[idx]
                #axs[x,y].scatter(xx,yy,c=z, s=0.1, cmap='viridis')# vmin=mindist, vmax=maxdist)
        
                if cl==None:
                    axs[x,y].hist2d(df_par[py],df_par[px],bins=15, norm = colors.LogNorm())
                else:
                	z=np.array(cl)
                	idx = z.argsort()
                	xx =df_par[px][idx]
                	yy = df_par[py][idx]
               		z = z[idx]
                	axs[x,y].scatter(xx,yy,c=z, s=0.1, cmap='viridis', norm = colors.LogNorm())# vmin=mindist, vmax=maxdist)

                axs[x,y].set_xlim(parlist[y]['lower_limit'],parlist[y]['upper_limit'])
                axs[x,y].set_ylim(parlist[x]['lower_limit'],parlist[x]['upper_limit'])

            
            niceaxis(axs,x,y,px,py,parlist,parlist,6)

    plt.subplots_adjust(wspace=None, hspace=None)

    return axs


#======================================================================================================================================================================================================================================================================

#plot for fitting here

#======================================================================================================================================================================================================================================================================



def compare_plot(p,filename,meq,datafile,modeltype,lw=0.4):
       # gmin,gmax,rmin,rmax=meq.Get_data4(datafile,p[0])
        gmax,gmin,rmax,rmin=meq.Get_data(datafile)
        A=gmin.index.values
        I=gmin.columns.values
        maxi= np.nanmax([ np.nanmax(rmax.to_numpy()),np.nanmax(gmax.to_numpy())])
        mini= np.nanmin([ np.nanmin(rmin.to_numpy()),np.nanmin(gmin.to_numpy())])

        fig, axs = plt.subplots(6, 2)

        for pi in p:
            ss=meq.findss(A,I,pi,modeltype)
            M=np.nanmax(ss[:,:,:,:],axis=2)
            m=np.nanmin(ss[:,:,:,:],axis=2)
          
            
            for ii,i in enumerate(I):
                axs[ii,0].plot(A,M[:,ii,0],'-',c=colorGREEN,linewidth=lw)
                axs[ii,1].plot(A,M[:,ii,1],'-',c=colorRED,linewidth=lw)    
                axs[ii,0].plot(A,m[:,ii,0],'--',c=colorCYAN,linewidth=lw)
                axs[ii,1].plot(A,m[:,ii,1],'--',c=colorPurple,linewidth=lw)
                              
                axs[ii,0].set_ylim(ymin=mini-0.15*mini,ymax=maxi+.15*maxi)
                axs[ii,1].set_ylim(ymin=mini-0.15*mini,ymax=maxi+.15*maxi)


               
        for ii,i in enumerate(I):

                axs[ii,0].plot(A,gmax.to_numpy()[:,ii],'ko', markersize=4.)
                axs[ii,0].plot(A,gmin.to_numpy()[:,ii],'ko', markersize=4., mfc='none')

                axs[ii,1].plot(A,rmax.to_numpy()[:,ii],'ko', markersize=4.)
                axs[ii,1].plot(A,rmin.to_numpy()[:,ii],'ko', markersize=4., mfc='none')





def bifu_plot_fit(par,A0,I0,axs,i,j,meq,model):
	#A0[0]=0
	colors = [colorBLUE, colorGREEN] # first color is black, last is red


	ss=meq.findss(A0,I0,par,model)
	sse=meq.getEigen(ss,A0,I0, par)
	

	for ai,a in enumerate(A0):
		for ipi,ip in enumerate(I0):

			for si in np.arange(0,ss.shape[2]):
				e=sse[ai,ipi,si,:]
				er=e.real
				ei=e.imag
				
				sss=ss[ai,ipi,si,:]


				if np.all(er!=0): 
					#print(sse)
					if np.all(er<0): #stable
						#if np.all(ei==0):
						#node
							m='o'
							axs[i+ipi,j].plot(np.log10(a),sss[0],m, c=colorGREEN ,markersize=1)
							axs[i+ipi,j+1].plot(np.log10(a),sss[2],m, c=colorRED, markersize=1.)
							#axs[i+ipi,j+2].plot(np.log10(a),sss[3],m, c= colorBLUE ,markersize=1.)


					elif np.any(er>0): #unstable
						#check complex conjugate
						if len(er) != len(set(er)):
							if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
								#oscilation
								m='o'
								axs[i+ipi,j+0].plot(np.log10(a),sss[0],m+'m',markersize=1.)
								axs[i+ipi,j+1].plot(np.log10(a),sss[2],m+'m',markersize=1.)
								#axs[i+ipi,j+2].plot(np.log10(a),sss[3],m+'m',markersize=1.)
							else:
								m='x'
								axs[i+ipi,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
								axs[i+ipi,j+1].plot(np.log10(a),sss[2],m+'k',markersize=1.)
								#axs[i+ipi,j+2].plot(np.log10(a),sss[3],m+'k',markersize=1.)
						else:
								m='x'
								axs[i+ipi,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
								axs[i+ipi,j+1].plot(np.log10(a),sss[2],m+'k',markersize=1.)
								#axs[i+ipi,j+2].plot(np.log10(a),sss[3],m+'k',markersize=1.)

	#axs[i,j].set_ylim(-0.1,1.1)
	#axs[i,j+1].set_ylim(-0.1,1.1)
	#axs[i,j+2].set_ylim(-0.1,1.1)

	return axs


def bifu_plot_fit_tgh(par,A0,I0,axs,i,j,meq,model):
	#A0[0]=0
	colors = [colorBLUE, colorGREEN] # first color is black, last is red


	ss=meq.findss(A0,I0,par,model)
	sse=meq.getEigen(ss,A0,I0, par)
	

	for ai,a in enumerate(A0):
		for ipi,ip in enumerate(I0):

			for si in np.arange(0,ss.shape[2]):
				e=sse[ai,ipi,si,:]
				er=e.real
				ei=e.imag
				
				sss=ss[ai,ipi,si,:]


				if np.all(er!=0): 
					#print(sse)
					if np.all(er<0): #stable
						#if np.all(ei==0):
						#node
							m='o'
							axs[i+ipi,j].plot(np.log10(a),sss[0],m, c=colorGREEN ,markersize=1)
							axs[i+ipi,j].plot(np.log10(a),sss[2],m, c=colorRED, markersize=1.)
							#axs[i+ipi,j].plot(np.log10(a),sss[3],m, c= colorBLUE ,markersize=1.)


					elif np.any(er>0): #unstable
						#check complex conjugate
						if len(er) != len(set(er)):
							if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
								#oscilation
								m='o'
								axs[i+ipi,j+0].plot(np.log10(a),sss[0],m+'m',markersize=1.)
								axs[i+ipi,j].plot(np.log10(a),sss[2],m+'m',markersize=1.)
								#axs[i+ipi,j].plot(np.log10(a),sss[3],m+'m',markersize=1.)
							else:
								m='x'
								axs[i+ipi,j].plot(np.log10(a),sss[0],m,c=colorGREEN,markersize=1.)
								axs[i+ipi,j].plot(np.log10(a),sss[2],m,c=colorRED,markersize=1.)
								#axs[i+ipi,j].plot(np.log10(a),sss[3],m+'k',markersize=1.)
						else:
								m='x'
								axs[i+ipi,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
								axs[i+ipi,j].plot(np.log10(a),sss[2],m+'k',markersize=1.)
								#axs[i+ipi,j].plot(np.log10(a),sss[3],m+'k',markersize=1.)

	#axs[i,j].set_ylim(-0.1,1.1)
	#axs[i,j+1].set_ylim(-0.1,1.1)
	#axs[i,j+2].set_ylim(-0.1,1.1)


	return axs


def bifu_plot_fit_tgh2(par,A0,I0,axs,i,j,meq,model):
	#A0[0]=0
	colors = [colorBLUE, colorGREEN] # first color is black, last is red


	ss=meq.findss(A0,I0,par,model)
	sse=meq.getEigen(ss,A0,I0, par)
	

	for ai,a in enumerate(A0):
		for ipi,ip in enumerate(I0):

			for si in np.arange(0,ss.shape[2]):
				e=sse[ai,ipi,si,:]
				er=e.real
				ei=e.imag
				
				sss=ss[ai,ipi,si,:]
				sss[0]=sss[0]-10**par['basal_green']


				if np.all(er!=0): 
					#print(sse)
					if np.all(er<0): #stable
						#if np.all(ei==0):
						#node
							m='o'
							axs[i+ai,j].plot(np.log10(ip),sss[0],m, c=colorGREEN ,markersize=1)
							axs[i+ai,j].plot(np.log10(ip),sss[2],m, c=colorRED, markersize=1.)
							#axs[i+ipi,j].plot(np.log10(a),sss[3],m, c= colorBLUE ,markersize=1.)


					elif np.any(er>0): #unstable
						#check complex conjugate
						if len(er) != len(set(er)):
							if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
								#oscilation
								m='o'
								axs[i+ai,j+0].plot(np.log10(ip),sss[0],m+'m',markersize=1.)
								axs[i+ai,j].plot(np.log10(ip),sss[2],m+'m',markersize=1.)
								#axs[i+ipi,j].plot(np.log10(a),sss[3],m+'m',markersize=1.)
							else:
								m='x'
								axs[i+ai,j].plot(np.log10(ip),sss[0],m,c=colorGREEN,markersize=1.)
								axs[i+ai,j].plot(np.log10(ip),sss[2],m,c=colorRED,markersize=1.)
								#axs[i+ipi,j].plot(np.log10(a),sss[3],m+'k',markersize=1.)
						else:
								m='x'
								axs[i+ai,j].plot(np.log10(ip),sss[0],m+'k',markersize=1.)
								axs[i+ai,j].plot(np.log10(ip),sss[2],m+'k',markersize=1.)
								#axs[i+ipi,j].plot(np.log10(a),sss[3],m+'k',markersize=1.)

	#axs[i,j].set_ylim(-0.1,1.1)
	#axs[i,j+1].set_ylim(-0.1,1.1)
	#axs[i,j+2].set_ylim(-0.1,1.1)

	return axs

def bifu_2Dplot_IPTG(A0,par,axs,i,j,IPTG,p,meq,modeltype,l=2):
	mono_matrix=np.zeros((len(A0),len(IPTG)))
	mono_green=np.zeros((len(A0),len(IPTG)))
	mono_red=np.zeros((len(A0),len(IPTG)))
	oscil_matrix=np.zeros((len(A0),len(IPTG)))

	maxvalues=1000#par['beta_green']
	minvalues=0.001#par['alpha_green']

	ss=meq.findss(A0,IPTG,par,modeltype)
	sse=meq.getEigen(ss,A0,IPTG,par,modeltype)

	oscil=np.apply_along_axis(isOscillation,3,sse)
	c_oscil=np.count_nonzero(oscil,axis=2)
	c=np.count_nonzero(~np.isnan(ss[:,:,:,0]),axis=2).astype(np.float)

	mono_matrix=c
	mono_green=ss[:,:,0,0]#-minvalues)/(maxvalues-minvalues)
	mono_red=ss[:,:,0,l]#-minvalues)/(maxvalues-minvalues)

	mono_green, mono_red = bestMatrix(mono_green,mono_red)
	oscil_matrix=c_oscil

	mono_matrix[mono_matrix==1]=np.nan
	#mono_green[mono_green<0.05*maxvalues]=np.nan
	#mono_red[mono_red<0.05*maxvalues]=np.nan
	oscil_matrix[oscil_matrix==0]=oscil_matrix[oscil_matrix==0]*np.nan

	
	axs[i,j].imshow(mono_red.T,aspect="auto", cmap="Reds", vmin=minvalues,vmax=maxvalues)
	axs[i,j].imshow(mono_green.T,aspect="auto", cmap="Greens", alpha=1,vmin=minvalues,vmax=maxvalues)
	axs[i,j].imshow(mono_matrix.T,aspect="auto", cmap="Blues", vmin=2,vmax=3)# vmin=1,vmax=5)
	#axs[i,j].imshow(oscil_matrix.T,aspect="auto", cmap="Greys", vmin=0,vmax=1)# vmin=1,vmax=5)
	a=np.argwhere(oscil_matrix.T==1)
	axs[i,j].plot(a[:,1],a[:,0],'mx',markersize=.5)

	indextick=[]
	for ind in np.arange(5):
		new_i = int(len(IPTG)/5)*ind
		indextick.append(new_i)

	#indextick=np.arange(0,len(k))
	Xlabel=np.linspace(min(IPTG),max(IPTG),5)
	Ylabel=np.logspace(min(A0),max(A0),5)

	axs[i,j].set_yticks(indextick)
	axs[i,j].set_xticks(indextick)
	axs[i,j].set_yticklabels(np.round(IPTG[indextick],2), fontsize=6 )
	axs[i,j].set_xticklabels(np.round(A0[indextick],3),rotation=90,fontsize=6 )

	'''
	axs[i,j].set_yticks(np.arange(len(k)))
	axs[i,j].set_xticks(np.arange(len(A0)))
	axs[i,j].set_yticklabels(np.round(k,2), fontsize=6)
	axs[i,j].set_xticklabels(np.round(A0,2),rotation=90,fontsize=6)
	'''

	return axs



def TI_plot_fit(AHL,IPTG,p0,meq):
	nm=100
	ee=meq.TuringInstability(AHL,IPTG,p0,nm,"TSXLT")
	fig, axs = plt.subplots(6,6,constrained_layout=True, figsize=(10,6))
	col= ['r','g','b','y','m']
	x=np.arange(nm)
	for ai,a in enumerate(AHL):
		for ii,i in enumerate(IPTG):

			for s in np.arange(5):


				axs[ai,ii].plot(x,ee.real[ai,ii,s,:,:],color=col[s])
			#axs[ai,ii].plot(x,ee.real[ai,ii,:,s,1].T)
			#axs[ai,ii].plot(x,ee.real[ai,ii,:,s,2].T)
			axs[ai,ii].set_yscale("symlog")
			axs[ai,ii].axhline(y=0,color='k')

			#axs[ai,ii].set_ylim(-0.00001,0.00001)
	return axs


'''



======================================================================================================================================================================================================================================================================

Script to handle multiple plot more nicely
need to fuse the two and improve them in one future

======================================================================================================================================================================================================================================================================

'''



def niceaxis(axs,x,y,px,py,kx,ky,size):
        #axs[x,y].set_ylim(-0.1,1.1)
        #if x==0:
        #    axs[x,y].set_title(py,fontsize=size)
        
        if y==0:
            axs[x,y].set_ylabel(px,fontsize=size,rotation=45,ha='right')
            axs[x,y].tick_params(axis='y', labelsize=size-2)

            if x!=len(kx)-1:
                    axs[x,y].set_xticks([])
                    axs[x,y].set_xticklabels([])
            #axs[x,y].set_yticks([])

        if x==len(kx)-1:
            axs[x,y].set_xlabel(py,fontsize=size,rotation=45)
            axs[x,y].tick_params(axis='x', labelsize=size-2)

            if y!=0:
                    axs[x,y].set_yticks([])
                    axs[x,y].set_yticklabels([])

        if y!=0 and x!=len(kx)-1:
            axs[x,y].set_xticks([])
            axs[x,y].set_yticks([])
            axs[x,y].set_xticklabels([])
            axs[x,y].set_yticklabels([])
        #if x==len(kx)-1:
        #   axs[x,y].set_xlabel('AHL')
        return axs


def niceaxis3(axs,x,y,p3,px,py,kx):
		#axs[x,y].set_ylim(-0.1,1.1)
		if x==0:
			axs[x,y].set_title(str(np.round(py,2)),fontsize=8)
		

		if y==0:
			axs[x,y].set_ylabel(p3,fontsize=8)
			if x!=len(kx)-1:
					axs[x,y].set_xticks([])
					axs[x,y].set_xticklabels([])
			#axs[x,y].set_yticks([])

		if x==len(kx)-1:
			if y!=0:
					axs[x,y].set_yticks([])
					axs[x,y].set_yticklabels([])

			
		if y==len(kx)-1:
			#axs[x,y].yaxis.tick_right()
			axs[x,y].set_yticks([])
			axs[x,y].set_yticklabels([])
			axs[x,y] = axs[x,y].twinx()
			axs[x,y].set_ylabel(str(np.round(px,2)),fontsize=8)
			axs[x,y].set_yticks([])
			axs[x,y].set_yticklabels([])
			if x!=len(kx)-1:
					axs[x,y].set_xticks([])
					axs[x,y].set_xticklabels([])


		if y!=0 and x!=len(kx)-1:
			axs[x,y].set_xticks([])
			axs[x,y].set_yticks([])
			axs[x,y].set_xticklabels([])
			axs[x,y].set_yticklabels([])
		#if x==len(kx)-1:
		#	axs[x,y].set_xlabel('AHL')
		return axs

def niceaxis2(axs,x,y,p3,px,py,kx,ky):
		#axs[x,y].set_ylim(-0.1,1.1)
		if x==0:
			axs[x,y].set_title(str(np.round(py,2)),fontsize=8)
		

		if y==0:
			axs[x,y].set_ylabel(str(np.round(px,3)),fontsize=8)
			if x!=len(kx)-1:
					axs[x,y].set_xticks([])
					axs[x,y].set_xticklabels([])
			#axs[x,y].set_yticks([])

		if x==len(kx)-1:
			if y!=0:
					axs[x,y].set_yticks([])
					axs[x,y].set_yticklabels([])

		'''		
		if y==len(ky):
			#axs[x,y].yaxis.tick_right()
			axs[x,y].set_yticks([])
			axs[x,y].set_yticklabels([])
			axs[x,y] = axs[x,y].twinx()
			axs[x,y].set_ylabel(str(np.round(px,2)),fontsize=8,fontname=font)
			axs[x,y].set_yticks([])
			axs[x,y].set_yticklabels([])
			if x!=len(kx)-1:
					axs[x,y].set_xticks([])
					axs[x,y].set_xticklabels([])
		'''

		if y!=0 and x!=len(kx)-1:
			axs[x,y].set_xticks([])
			axs[x,y].set_yticks([])
			axs[x,y].set_xticklabels([])
			axs[x,y].set_yticklabels([])
		#if x==len(kx)-1:
		#	axs[x,y].set_xlabel('AHL')
		return axs




'''
==================================================================================================================================================================================================================================================================



other


======================================================================================================================================================================================================================================================================
'''


def load(number, filename,parlist):
    namelist=[]
    for i,par in enumerate(parlist):
        namelist.append(parlist[i]['name'])
    
    path = filename+'/smc/pars_' + number + '.out'
    dist_path = filename+'/smc/distances_' + number + '.out'

    raw_output= np.loadtxt(path)
    dist_output= np.loadtxt(dist_path)
    df = pd.DataFrame(raw_output, columns = namelist)
    df['dist']=dist_output
    idx=np.argsort(df['dist'])
    df=df.sort_values('dist',ascending=False)
    p=[]

    for i in idx:
        p0=df.loc[i].tolist()
        p.append(pars_to_dict(p0,parlist))

    
    return p, df 


def pars_to_dict(pars,parlist):
### This function is not necessary, but it makes the code a bit easier to read,
### it transforms an array of pars e.g. p[0],p[1],p[2] into a
### named dictionary e.g. p['k0'],p['B'],p['n'],p['x0']
### so it is easier to follow the parameters in the code
    dict_pars = {}
    for ipar,par in enumerate(parlist):
        dict_pars[par['name']] = pars[ipar] 
    return dict_pars
'''
=========================================================================================================================================================================================================================================
old script in case
=======================================================================================================================================================================================================================================
'''




'''
def bifu_plot2(par):
	A0=np.logspace(-3,0,width)
	A0[0]=0


	ss=findss(A0,[0],par)


	#ss=np.sort(ss,axis=2)
	fig, axs = plt.subplots(3)
	for mi,m in enumerate(['-','--','-','o','o']):
                    axs[0].plot(ss[:,0,mi,0],m+'g',markersize=1.)
                    axs[1].plot(ss[:,0,mi,1],m+'r',markersize=1.)
                    axs[2].plot(ss[:,0,mi,2],m+'b',markersize=1.)

	plt.show()

def bifu_plot(A0,par):
	A0=np.logspace(-3,0,width)
	#A0[0]=0
	ss=findss(A0,[0],par)
	sse=getEigen(ss,A0,par)
	I=0

	fig, axs = plt.subplots(3)

	for ai,a in enumerate(A0):
		for si in np.arange(0,ss.shape[2]):
			e=sse[ai,I,si,:]
			er=e.real
			ei=e.imag
			

			sss=ss[ai,I,si,:]
			if np.all(er!=0): 
				#print(sse)
				if np.all(er<0): #stable
					#if np.all(ei==0):
					#node
						m='o'
						axs[0].plot(np.log10(a),sss[0],m+'g',markersize=1)
						axs[1].plot(np.log10(a),sss[1],m+'r',markersize=1.)
						axs[2].plot(np.log10(a),sss[2],m+'b',markersize=1.)


				elif np.any(er>0): #unstable
					#check complex conjugate
					if len(er) != len(set(er)):
						if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
							#oscilation
							m='o'
							axs[0].plot(np.log10(a),sss[0],m+'m',markersize=1.)
							axs[1].plot(np.log10(a),sss[1],m+'m',markersize=1.)
							axs[2].plot(np.log10(a),sss[2],m+'m',markersize=1.)
						else:
							m='x'
							axs[0].plot(np.log10(a),sss[0],m+'k',markersize=1.)
							axs[1].plot(np.log10(a),sss[1],m+'k',markersize=1.)
							axs[2].plot(np.log10(a),sss[2],m+'k',markersize=1.)
					else:
							m='x'
							axs[0].plot(np.log10(a),sss[0],m+'k',markersize=1.)
							axs[1].plot(np.log10(a),sss[1],m+'k',markersize=1.)
							axs[2].plot(np.log10(a),sss[2],m+'k',markersize=1.)


	plt.show()


def plot1D(par,model = "TSXLT", endpoint=False):


	n=np.arange(0,10)
	tt=int(time/dt/4)

	if endpoint==True:
		n=[1]
	fig, axs = plt.subplots(len(n),6,sharey=False)


	for v in np.arange(0,6):
		A0=np.logspace(-3,0,width)*0
		#A0=np.ones(width)*0.7


		if v==0:
			G=np.random.normal(loc=0.8, scale=0.2, size=(width))
			R=np.random.normal(loc=0.8, scale=0.2, size=(width))*0
			A=np.random.normal(loc=0.8, scale=0.2, size=(width))

		elif v==1:
			G=np.random.normal(loc=0.8, scale=0.2, size=(width))*0
			R=np.random.normal(loc=0.8, scale=0.2, size=(width))
			A=np.random.normal(loc=0.8, scale=0.2, size=(width))
		elif v==2:
			G=np.random.normal(loc=0.8, scale=0.2, size=(width))
			R=np.random.normal(loc=0.8, scale=0.2, size=(width))
			A=np.random.normal(loc=0.8, scale=0.2, size=(width))
		elif v==3:
			G=np.random.normal(loc=0.1, scale=0.1, size=(width))
			R=np.random.normal(loc=0.1, scale=0.1, size=(width))*0
			A=np.random.normal(loc=0.1, scale=0.1
				, size=(width))
		elif v==4:
			G=np.random.normal(loc=0.1, scale=0.1, size=(width))*0
			R=np.random.normal(loc=0.1, scale=0.1, size=(width))
			A=np.random.normal(loc=0.1, scale=0.1, size=(width))
		elif v==5:
			G=np.random.normal(loc=0.1, scale=0.1, size=(width))
			R=np.random.normal(loc=0.1, scale=0.1, size=(width))
			A=np.random.normal(loc=0.1, scale=0.1, size=(width))

		av=simulation(G,R,A,A0,0,par,dt,model,tt=time)
		if endpoint==False:
			for i in n:
				tt=int(time/dt/len(n))*i
				index=v
				axs[i,index].plot(av[tt,0],'g')
				axs[i,index].plot(av[tt,1],'r')
				#axs[i,index+2].plot(av[tt,2],'b')

				#axs[i,v].set_ylim(0,0.02)
		else:
				axs[v].imshow(av[-1,0])
				#axs[v].imshow(av[-1,1],'--r')


	#plt.plot(at[round(tt/2)])
	#plt.plot(at[-1])
	plt.show()



def plot2D(par,model = "TSXLT", endpoint=False):


	n=np.arange(0,10)
	tt=int(time/dt/4)

	if endpoint==True:
		n=[1]
	fig, axs = plt.subplots(len(n),6*3,sharey=False)


	for v in np.arange(0,6):
		A0=np.logspace(-3,0,width,width)*0
		#A0=np.ones(width)*0.7


		if v==0:
			G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
			R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
			A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))

		elif v==1:
			G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
			R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
			A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
		elif v==2:
			G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
			R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
			A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
		elif v==3:
			G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
			A=np.random.normal(loc=0.1, scale=0.1
				, size=(width,width))
		elif v==4:
			G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
			R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
		elif v==5:
			G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))

		av=simulation(G,R,A,A0,0,par,dt,model,tt=time)
		if endpoint==False:
			for i in n:
				tt=int(time/dt/len(n))*i
				index=3*v
				axs[i,index].imshow(av[tt,0],'Greens')
				axs[i,index+1].imshow(av[tt,1],'Reds')
				axs[i,index+2].imshow(av[tt,2],'Blues')

				#axs[i,v].set_ylim(0,0.02)
		else:
				axs[v].imshow(av[-1,0])
				#axs[v].imshow(av[-1,1],'--r')


	#plt.plot(at[round(tt/2)])
	#plt.plot(at[-1])
	plt.show()

def plot_eigen(par):
	A=np.zeros(1)
	e=get_turingInstability(A,par,model="TSXLT")
	er=e.real
	ei=e.imag
	plt.plot(er[0,0,0,:,:],'g')
	plt.plot(er[0,0,1,:,:],'r')
	plt.plot(er[0,0,2,:,:],'b')
	plt.plot(er[0,0,3,:,:],'y')
	plt.plot(er[0,0,4,:,:],'c')
	plt.yscale('symlog')
	plt.ylim(-100,1)
	plt.axhline(y = 0, linestyle='-',c='black')
	plt.show()

def plot_all(par,model,k=[0],parm="null",i=""):
	if len(k)==1:
		fig, axs = plt.subplots(1,10,sharey=False,figsize=(40,1))
		A0=np.logspace(-3,0,width)
		I=0

		axs[0].axis([0, 10, 0, 10])
		axs[0].text(1, 5, k, fontsize=5)

		#bifu plot

		ss=findss(A0,[0],par)
		sse=getEigen(ss,A0,par)

		for ai,a in enumerate(A0):
			for si in np.arange(0,ss.shape[2]):
				e=sse[ai,I,si,:]
				er=e.real
				ei=e.imag
				

				sss=ss[ai,I,si,:]
				if np.all(er!=0): 
					#print(sse)
					if np.all(er<0): #stable
						#if np.all(ei==0):
						#node
							m='o'
							axs[1].plot(np.log10(a),sss[0],m+'g',markersize=1)
							axs[2].plot(np.log10(a),sss[1],m+'r',markersize=1.)
							#axs[0].plot(np.log10(a),sss[2],m+'b',markersize=1.)


					elif np.any(er>0): #unstable
						#check complex conjugate
						if len(er) != len(set(er)):
							if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
								#oscilation
								m='o'
								axs[1].plot(np.log10(a),sss[0],m+'m',markersize=1.)
								axs[2].plot(np.log10(a),sss[1],m+'m',markersize=1.)
								#axs[0].plot(np.log10(a),sss[2],m+'m',markersize=1.)
							else:
								m='x'
								axs[1].plot(np.log10(a),sss[0],m+'k',markersize=1.)
								axs[2].plot(np.log10(a),sss[1],m+'k',markersize=1.)
								#axs[0].plot(np.log10(a),sss[2],m+'k',markersize=1.)
						else:
								m='x'
								axs[1].plot(np.log10(a),sss[0],m+'k',markersize=1.)
								axs[2].plot(np.log10(a),sss[1],m+'k',markersize=1.)
								#axs[0].plot(np.log10(a),sss[2],m+'k',markersize=1.)
		#1D plot


		for v in np.arange(0,6):
			A0=np.logspace(-3,0,width,width)*0

			if v==0:
				G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
				A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))

			elif v==1:
				G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
				R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
			elif v==2:
				G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
			elif v==3:
				G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
				A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			elif v==4:
				G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
				R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			elif v==5:
				G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))

			av=simulation(G,R,A,A0,0,par,dt,model,tt=time)

			axs[v+3].imshow(av[-1,0],cmap='Greens')
			#axs[v+3].imshow(av[-1,1],cmap='Reds')

		#eigen
		A=np.zeros(1)
		e=get_turingInstability(A,par)
		er=e.real
		ei=e.imag
		npl = 9
		axs[npl].plot(er[0,0,0,:,:],'g')
		axs[npl].plot(er[0,0,1,:,:],'r')
		axs[npl].plot(er[0,0,2,:,:],'b')
		axs[npl].plot(er[0,0,3,:,:],'y')
		axs[npl].plot(er[0,0,4,:,:],'c')
		axs[npl].set_yscale('symlog')
		axs[npl].set_ylim(-100,1)
		axs[npl].axhline(y = 0, linestyle='-',c='black')

	else:
		fig, axs = plt.subplots(len(k),10,sharey=False,figsize=(10,len(k)))
		I=0

		#partxt = json.dumps(par)
		#partxt=str(par)
		#partxt=partxt.replace(",", "\n")
		#print(partxt)

		for ki,kii in enumerate(k):
			A0=np.logspace(-3,0,width)

			par[parm]=kii
			txt=parm + "\n" + str(kii)
			axs[ki,0].axis([0, 10, 0, 10])
			axs[ki,0].text(0.1, 5, txt, fontsize=8)
			#bifu plot
			ss=findss(A0,[0],par,model)
			sse=getEigen(ss,A0,par,model)

			for ai,a in enumerate(A0):
				for si in np.arange(0,ss.shape[2]):
					e=sse[ai,I,si,:]
					er=e.real
					ei=e.imag

					sss=ss[ai,I,si,:]
					if np.all(er!=0): 
						#print(sse)
						if np.all(er<0): #stable
							#if np.all(ei==0):
							#node
								m='o'
								axs[ki,1].plot(np.log10(a),sss[0],m+'g',markersize=1)
								axs[ki,2].plot(np.log10(a),sss[1],m+'r',markersize=1.)
								#axs[0].plot(np.log10(a),sss[2],m+'b',markersize=1.)


						elif np.any(er>0): #unstable
							#check complex conjugate
							if len(er) != len(set(er)):
								if np.any(ei==0) and len(ei) != len(set(np.abs(ei))):
									#oscilation
									m='o'
									axs[ki,1].plot(np.log10(a),sss[0],m+'m',markersize=1.)
									axs[ki,2].plot(np.log10(a),sss[1],m+'m',markersize=1.)
									#axs[0].plot(np.log10(a),sss[2],m+'m',markersize=1.)
								else:
									m='x'
									axs[ki,1].plot(np.log10(a),sss[0],m+'k',markersize=1.)
									axs[ki,2].plot(np.log10(a),sss[1],m+'k',markersize=1.)
									#axs[0].plot(np.log10(a),sss[2],m+'k',markersize=1.)
							else:
									m='x'
									axs[ki,1].plot(np.log10(a),sss[0],m+'k',markersize=1.)
									axs[ki,2].plot(np.log10(a),sss[1],m+'k',markersize=1.)
									#axs[0].plot(np.log10(a),sss[2],m+'k',markersize=1.)
			
			axs[ki,1].set_ylim(-0.1,1.1)
			axs[ki,2].set_ylim(-0.1,1.1)



			#1D plot


			for v in np.arange(0,6):
				A0=np.logspace(-3,0,width,width)*0

				if v==0:
					G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
					A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))

				elif v==1:
					G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
					R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				elif v==2:
					G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				elif v==3:
					G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
					A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				elif v==4:
					G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
					R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				elif v==5:
					G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))

				av=simulation(G,R,A,A0,0,par,dt,model,tt=time)

				axs[ki,v+3].imshow(av[-1,0])
				#axs[ki,v+3].plot(av[-1,1],'--r')
				#axs[ki,v+3].set_ylim(-0.1,1.1)


			#eigen
			A=np.zeros(1)
			e=get_turingInstability(A,par)
			er=e.real
			ei=e.imag
			npl = 9
			axs[ki,npl].plot(er[0,0,0,:,:],'g')
			axs[ki,npl].plot(er[0,0,1,:,:],'r')
			axs[ki,npl].plot(er[0,0,2,:,:],'b')
			axs[ki,npl].plot(er[0,0,3,:,:],'y')
			axs[ki,npl].plot(er[0,0,4,:,:],'c')
			#axs[ki,npl].set_yscale('symlog')
			axs[ki,npl].set_ylim(-0.5,0.5)
			axs[ki,npl].axhline(y = 0, linestyle='-',c='black')
	plt.savefig(parm+'_'+i+'.png')
	plt.show()



def screen_2par(par,par1,par2,k1,k2):
	size=len(k1)

	I=0
	A0=np.logspace(-2,0,20)



	fig, axs = plt.subplots(size,size,sharey='none') #figsize=(12,8)
	for i,ii in enumerate(k1):
		for j,jj in enumerate(k2):
			par[par1]=ii#ii
			par[par2]= jj#ii#ii

			ss=findss(A0,[0],par)
			sse=getEigen(ss,A0,par)

			for ai,a in enumerate(A0):
					for si in np.arange(0,ss.shape[2]):
						sss=ss[ai,I,si,:]
						e=sse[ai,I,si,:]

						if np.all(e!=0): 
							#print(sse)
							if np.all(e<0):
								m='o'
								axs[i,j].plot(np.log10(a),sss[0],m+'g',markersize=1)
								#axs[i,j].plot(np.log10(a),sss[1],m+'r',markersize=1.)
								#axs[2].plot(np.log10(a),sss[2],m+'b',markersize=1.)

							elif np.any(e>0):
									if len(e) != len(set(e)):
										m='o'
										axs[i,j].plot(np.log10(a),sss[0],m+'m',markersize=1.)
										#axs[i,j].plot(np.log10(a),sss[1],m+'m',markersize=1.)
										#axs[2].plot(np.log10(a),sss[2],m+'m',markersize=1.)
									else:
										m='x'
										axs[i,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
										#axs[i,j].plot(np.log10(a),sss[1],m+'k',markersize=1.)
										#axs[2].plot(np.log10(a),sss[2],m+'k',markersize=1.)
								
							else:
									m='x'
									axs[i,j].plot(np.log10(a),sss[0],m+'k',markersize=1.)
									#axs[i,j].plot(np.log10(a),sss[1],m+'k',markersize=1.)
									#axs[2].plot(np.log10(a),sss[2],m+'k',markersize=1.)

			axs[i,j].axes.xaxis.set_ticks([])
			axs[i,j].axes.yaxis.set_ticks([])
			axs[i,j].set_ylim(-0.1,1.1)
	plt.show()

'''