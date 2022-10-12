from model import *
from diffusion import *
from grow import *


import pandas as pd
import sys

from scipy.stats import gaussian_kde
import numpy as np

def simulation(G,R,A,C,D,A0,IPTG,par,width,dx,time,dt,model = "TSXLT",tt=1,noise=0,stationary=True,store=True,meq=None):
	ti=1
	#av=np.zeros((int(tt/dt),3,len(G)))

	Gi=np.copy(G)
	Ri=np.copy(R)
	Rfi=np.copy(R)

	Ai=np.copy(A)
	Ci=np.copy(C)
	Di=np.copy(D)

	center = np.where(C>0)
	gcell=[[],[]]

	if store:
		av=np.zeros((int(tt/dt),6,len(G),len(G)))
		av[0,0]=Gi
		av[0,2]=Rfi
		av[0,1]=Ri
		av[0,3]=Ai
		av[0,4]=Ci
		av[0,5]=Di
	else:
		av=np.zeros((1,6,len(G),len(G)))


	I=IPTG

	while (ti < int(tt/dt)):


			#1. make cell growth

			if stationary:
					stat_effect = calculateStationary(Di,par_growth) 
			else:
					stat_effect = 0

			d =  densityGrowth_model(Di,par_growth)


			#2. extend colonies size if less than 6 colonies on simulation

			if len(center[0])<6:
					gcell =  horizontalGrowth_model(Ci,center,ti*dt,dx,par_growth)
			else:
				gcell=[]

			#3. calculated state/fluo/diffusion
			if model== "TSXLT":

				#3.1 diffusion
				'''
				#change suze of matrix for diffusion
				#Ai_reduce=decrease_matrix(Ai,dx,dx_diffusion)
				#adif_small = diffusion2D(Ai_reduce,par['D_ahl'],dx_diffusion,dt)
				#adif = increase_matrix(adif_small,dx,dx_diffusion)
				'''
				adif=diffusion2D(Ai,par['D_ahl'],dx,dt)

				#3.2 dX/dt
				g,r,a =  model_TSXLT(Gi,Ri,Ai,A0,I,par)
				Gi += dt*g*(Ci - stat_effect)
				Ri += dt*r*(Ci - stat_effect)
				Ai += dt*a*(Ci - stat_effect) + adif

			elif model== "TSXLT_fit":

				#3.1 diffusion
				'''
				#change suze of matrix for diffusion
				#Ai_reduce=decrease_matrix(Ai,dx,dx_diffusion)
				#adif_small = diffusion2D(Ai_reduce,par['D_ahl'],dx_diffusion,dt)
				#adif = increase_matrix(adif_small,dx,dx_diffusion)
				'''
				adif=diffusion2D(Ai,par['D_ahl'],dx,dt)

				#3.2 dX/dt

				g,r,rf,a =  meq.model_TSXLT(Gi,Ri,Rfi,Ai,A0,I,par)
				#noise
			


				Gi += dt*g*(Ci - stat_effect) 
				Ri += dt*r*(Ci - stat_effect)  
				Rfi += dt*rf*(Ci - stat_effect)
				#Ai += dt*a*(Ci - stat_effect)*Di/par_growth['max_density'] + adif
				Ai += dt*a*(Ci - stat_effect) + adif
				
				if noise>0:
					nnn=noise
					gn= (np.random.normal(loc=0, scale=nnn, size=(width,width)))*(Ci - stat_effect)
					rn= (np.random.normal(loc=0, scale=nnn, size=(width,width)))*(Ci - stat_effect)
					#an= (np.random.normal(loc=0, scale=nnn, size=(width,width)))*(Ci - stat_effect)

					Gi+=gn
					Ri+=rn
				#Ai+=an

				Gi[Gi<0]=0
				Ri[Ri<0]=0



			elif model== "Turing":
				g,r = model_turing(Gi,Ri,par)
				Gi += dt*g + diffusion2D(Gi,par['D_green'],dx,dt) #+gh
				Ri += dt*r + diffusion2D(Ri,par['D_red'],dx,dt) #+rh
				av[ti,0]=Gi 
				av[ti,1]=Ri

			elif model== "Repressilator":

				g,r,b = model_Repressilator(Gi,Ri,Ai,par)

				Gi += dt*g*(Ci - stat_effect)
				Ri += dt*r*(Ci - stat_effect)
				Ai += dt*b*(Ci - stat_effect)

				#gh,rh,bh =heredity(Gi,Ri,Ai,Ci,gcell,width)
				#Ci[gcell[:,0,0],gcell[:,0,1]]=1
				#Di[gcell[:,0,0],gcell[:,0,1]]=1
				#Gi += gh
				#Ri += rh
				#Ai += bh
				#av[ti,0]=Gi 
				#av[ti,1]=Ri
				#av[ti,3]=Ai 
				#av[ti,4]=Ci
				#av[ti,5]=Di

			#4. add new cells
			if len(gcell)>0:
				gh,rh,ah =heredity(Gi,Ri,Ai,Ci,gcell,width)	
				Gi += gh
				Ri += rh
				Ai += ah
				Ci[gcell[:,0,0],gcell[:,0,1]]=1
				Di[gcell[:,0,0],gcell[:,0,1]]=1

			#5.store values
			Di =  Di + dt*d

			if store:
				av[ti,0]=Gi 
				av[ti,1]=Ri 
				av[ti,2]=Rfi 
				av[ti,3]=Ai 
				av[ti,4]=Ci
				av[ti,5]=Di 

			ti += 1
	if store == False:
		av[0,0]=Gi 
		av[0,1]=Ri 
		av[0,2]=Rfi 
		av[0,3]=Ai 
		av[0,4]=Ci
		av[0,5]=Di 
	return av


def init_grid(values,noise,width,Tmatrix):
	G=np.random.normal(loc=values[0], scale=noise[0], size=(width,width))*Tmatrix[0]
	R=np.random.normal(loc=values[1], scale=noise[1], size=(width,width))*Tmatrix[1]
	A=np.random.normal(loc=values[2], scale=noise[2], size=(width,width))*Tmatrix[2]
	D=np.random.normal(loc=values[3], scale=noise[3], size=(width,width))
	C=np.ones((width,width))
	return G,R,A,C,D


def init_colony(values,width):
	G=np.zeros((width,width))
	R=np.zeros((width,width))
	A=np.zeros((width,width))
	D=np.zeros((width,width))
	C=np.zeros((width,width))
	C[int(width/2),int(width/2)]=1
	D[int(width/2),int(width/2)]=1
	G[int(width/2),int(width/2)]=values[0]
	R[int(width/2),int(width/2)]=values[1]
	A[int(width/2),int(width/2)]=values[2]

	return G,R,A,C,D





def getMAP(points):
	kernel_estimate = gaussian_kde(points) # generate a continuous kernel density from the data
	densities = kernel_estimate(points) # evaluate the data on the kernel
	return points.T[np.argmax(densities)] # return the parameter set that is more dense




def pars_to_dict(pars,parlist):
### This function is not necessary, but it makes the code a bit easier to read,
### it transforms an array of pars e.g. p[0],p[1],p[2] into a
### named dictionary e.g. p['k0'],p['B'],p['n'],p['x0']
### so it is easier to follow the parameters in the code
    dict_pars = {}
    for ipar,par in enumerate(parlist):
        dict_pars[par['name']] = pars[ipar] 
    return dict_pars
##plotting part