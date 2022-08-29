
'''
Here are the script to run stuff

written by Içvara
'''
from model import *
from diffusion import *
from grow import *
from  easyplot import *
from TuringInstability import*
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import sys


#check dx???


def simulation(G,R,A,C,D,A0,IPTG,par,width,dx,time,dt,model = "TSXLT",tt=1,stationary=True):
	ti=1
	av=np.zeros((int(tt/dt),6,len(G),len(G)))
	#av=np.zeros((int(tt/dt),3,len(G)))

	Gi=np.copy(G)
	Ri=np.copy(R)
	Rfi=np.copy(R)

	Ai=np.copy(A)
	Ci=np.copy(C)
	Di=np.copy(D)

	center = np.where(C>0)
	gcell=[[],[]]

	av[0,0]=Gi
	av[0,1]=Rfi
	av[0,2]=Ri
	av[0,3]=Ai
	av[0,4]=Ci
	av[0,5]=Di

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
				Gi += dt*g*(Ci - stat_effect)
				Ri += dt*r*(Ci - stat_effect)
				Rfi += dt*rf*(Ci - stat_effect)
				Ai += dt*a*(Ci - stat_effect)*Di/par_growth['max_density'] + adif
				#Ai += dt*a*(Ci - stat_effect) + adif


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

				gh,rh,bh =heredity(Gi,Ri,Ai,Ci,gcell,width)
				Ci[gcell[:,0,0],gcell[:,0,1]]=1
				Di[gcell[:,0,0],gcell[:,0,1]]=1
				Gi += gh
				Ri += rh
				Ai += bh
				av[ti,0]=Gi 
				av[ti,1]=Ri
				av[ti,2]=Ai 
				av[ti,3]=Ci
				av[ti,4]=Di

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

			av[ti,0]=Gi 
			av[ti,1]=Rfi 
			av[ti,2]=Ri 
			av[ti,3]=Ai 
			av[ti,4]=Ci
			av[ti,5]=Di 
 


			ti += 1

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

'''
================================================================================================================00

par setup HERE

======================================================================================================================

'''


par0={
	'beta_green':1,
	'alpha_green':.001,#0.001,

	'K_ahl_green':0,
	'n_ahl_green':2,
	'delta_green':1,

	'K_GREEN':1,
	'n_GREEN':2,

	'beta_red':1,
	'alpha_red':0.001,

	'K_ahl_red':0,
	'n_ahl_red':2,
	'delta_red':1.,

	'K_RED':1,
	'n_RED':2,

	'beta_ahl':0,
	'K_ahl':0,
	'n_ahl':2,
	'delta_ahl':1,

    'D_ahl':1,#0.1,

    'K_IPTG':1
    }

size=2 #mm 10-1
dx=0.05
width=int(size/dx)

#dt= 0.9 * (dx ** 2) / (2 * D)

#ar0['D_ahl']=.01
dx_diffusion=0.05
time=16
dt=0.005



'''
================================================

Bifurcation plot and anlysis

===============================================
'''

'''
nPar1=100

A0=np.logspace(-3,0,nPar1)
k = np.linspace(0.,2,nPar1) #for K_RED/GREEN
p3='K_GREEN'



par=par0
#par['n_ahl_red']=0.5


#par['K_RED']=1.67

fig, axs = plt.subplots(3,2,constrained_layout=True)#,figsize=(1,1))
par=par0
par['beta_red']=1
par['alpha_red']=0.001
par['K_ahl_red']=1
par['beta_green']=0
par['alpha_green']=1
monostable_plot(A0,par,axs,0,0,k,p3)
axs[0,0].set_title("TSL")

par=par0
par['beta_red']=1
par['alpha_red']=0.001
par['K_ahl_green']=1
par['K_ahl_green']=1.5
par['K_ahl_red']=1
par['beta_green']=1
par['alpha_green']=0.001
monostable_plot(A0,par,axs,1,0,k,p3)
axs[1,0].set_title("TSLT")

par=par0
par['beta_red']=1
par['alpha_red']=0.001
par['K_ahl_green']=1
par['K_ahl_green']=1.5
par['K_ahl_red']=1
par['K_ahl_red']=1.25
par['beta_green']=1
par['alpha_green']=0.001
par['beta_ahl']=1
par['K_ahl']=0
monostable_plot(A0,par,axs,1,1,k,p3)
axs[1,1].set_title("TSXLT")

par=par0
par['beta_red']=1
par['alpha_red']=0.001
par['K_ahl_green']=1
par['K_ahl_green']=1.5
par['K_ahl_red']=1
par['beta_green']=0
par['alpha_green']=1
par['beta_ahl']=1
par['K_ahl']=0
monostable_plot(A0,par,axs,0,1,k,p3)
axs[0,1].set_title("TSXL")

par=par0
par['beta_red']=0
par['alpha_red']=1
par['K_ahl_green']=1
par['K_ahl_red']=1
par['beta_green']=1
par['alpha_green']=0.001
par['beta_ahl']=1
par['K_ahl']=0
monostable_plot(A0,par,axs,2,1,k,p3)
axs[2,1].set_title("TSXT")


par=par_TSLT_mushroom

monostable_plot(A0,par,axs,2,0,k,p3)
axs[2,0].set_title("TSLT_mushroom")

plt.savefig("Bifuplot_tilted.pdf",dpi=300)
plt.show()

'''


'''
old one
par=par_TSL
monostable_plot(A0,par,axs,0,0,k,p3)

par=par_TSXL
monostable_plot(A0,par,axs,1,0,k,p3)

par=par_TSLT
monostable_plot(A0,par,axs,0,1,k,p3)

par=par_TSXLT
monostable_plot(A0,par,axs,1,1,k,p3)

par=par_TSLT_mushroom
monostable_plot(A0,par,axs,0,2,k,p3)

par=par_TSXLT_mushroom
monostable_plot(A0,par,axs,1,2,k,p3)

par=par_TSXLT_mushroom2
monostable_plot(A0,par,axs,1,3,k,p3)


plt.savefig("Bifuplot.pdf",dpi=300)
plt.show()

stop
'''


'''
nPar1=100
A0=np.logspace(-3,0,nPar1)

par=par0
par['beta_red']=1
par['alpha_red']=0.001
par['K_ahl_green']=1
#par['K_ahl_green']=1.5
par['K_ahl_red']=1
#par['K_ahl_red']=1.25
par['beta_green']=1
par['alpha_green']=0.001
par['beta_ahl']=1
par['K_ahl']=0

#par=par_TSLT_mushroom

nPar=9
pKG=np.linspace(0,2,nPar)
#pKG=np.linspace(0.5,1.,nPar)


fig, axs = plt.subplots(nPar,3,figsize=(4,6),constrained_layout=True)

for pi,pG in enumerate(pKG):
	par['K_GREEN']=pG
	bifu_plot(par,A0,axs,pi,0)
	niceaxis(axs,pi,0,pG,'x',pG,pKG)

#plt.tight_layout()
plt.savefig("bifTSXLT1.pdf",dpi=300, transparent=False)
plt.show()


#par['beta_ahl']=1

stop
'''


#===============================================================================================================================================================


#screen of par space here


#======================================================================================================================================================================

'''
par=par0
par=par0
par['beta_red']=1
par['alpha_red']=0.5

par['K_ahl_green']=1
par['K_ahl_green']=1.75
par['K_ahl_red']=1
par['beta_green']=2
par['alpha_green']=0.001

name="TSXL_mush"


nPar1=30

A0=np.logspace(-3,0,nPar1)
k = np.linspace(0.,2,nPar1) #for K_RED/GREEN
p3='K_GREEN'

nPar=9

#try ratio n and k ahl
ky=np.linspace(1,3,nPar)
p1='K_ahl_green'
#delta red=1.67
#more n ar and lower ahl
kx=np.linspace(0,2,nPar)
p2='beta_red'

fig, axs = plt.subplots(len(kx),len(ky),constrained_layout=True)#,figsize=(1,1))

for x,px in enumerate(kx):
	for y,py in enumerate(ky):
		#for i,pk in enumerate(k):
		par[p1]=py
		par[p2]=px
			#par['K_GREEN']=pk
			#bifu_plot_par(par,A0,axs,x,y,i,nPar)
		monostable_plot(A0,par,axs,x,y,k,p3)
		niceaxis(axs,x,y,p3,px,py,kx)
		#axs[x,y].set_aspect('equal')


fig.supxlabel('AHL')
fig.supylabel(p2)
fig.suptitle(p1)

#plt.tight_layout()
#plt.subplots_adjust(wspace=0.1, hspace=0.1)

plt.savefig(name+"_"+p1+"_"+p2+"_"+p3+".pdf",dpi=300)
#plt.savefig("TEST.png",dpi=300, transparent=True)

plt.show()


stop
'''

#bifu plot
'''
fig, axs = plt.subplots(len(kx),len(ky))

for x,px in enumerate(kx):
	for y,py in enumerate(ky):
		for i,pk in enumerate(k):
				par[p1]=py
				par[p2]=px
				par[p3]=pk
				bifu_plot_par(par,A0,axs,x,y,i,len(k))
		#monostable_plot(A0,par,axs,x,y,k,'K_ahl_green')

		axs[x,y].set_ylim(-0.1,1.1)
		if x==0:
			axs[x,y].set_title(str(py))
		if y==0:
			axs[x,y].set_ylabel(str(px))
		if x==len(kx)-1:
			axs[x,y].set_xlabel('AHL')

plt.show()
'''


'''
================================================================================================================00

2D simulation HERE

======================================================================================================================
'''

'''
size=10 #mm 10-1
dx=0.1
width=int(size/dx)

 


dx_diffusion= dx
time=10
#dt=0.05

par=par0
par['beta_red']=1
par['alpha_red']=0.001
par['K_ahl_green']=1
par['K_ahl_green']=1.5
par['K_ahl_red']=1
par['beta_green']=0
par['alpha_green']=1
par['beta_ahl']=1
par['K_ahl']=0

#par=par_TSLT_mushroom

par['D_ahl']=.000001 #mm/h
par['D_ahl']=size/10/time #mm/h


Tmatrx=[1,1,0]
Vmatrx=[0,1.,1,1]
Nmatrx=[0.2,0.2,0.2,0.]

G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
#G,R,A,C,D = init_colony(Vmatrx,width)


name='TSXL_R_lawn2'



p1='K_GREEN'
kx = np.linspace(0.8,1.5,6)  #0,2
#kx = np.linspace(0.,0.5,6)  #0,2

#kx=np.array([0,0.8])


par['K_GREEN']=0.64

#A0
ky = np.logspace(-3,0.,10) #-3,0

dt =np.min( ((0.15 *dx**2 / par['D_ahl']) , 0.1))
print(dt)


#A0=np.ones((width,width))*A0_bifu[:,np.newaxis]


A0_bifu=np.logspace(-3,0.,width)

fig, axs = plt.subplots(len(kx),len(ky)+1,constrained_layout=True)

for x,kpx in enumerate(kx):
	par[p1]=kpx
	#Nmatrx=[kpx,kpx,0.,0.]
	G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
	bifu_plot_tgh(par,A0_bifu,axs,x,0)
	niceaxis2(axs,x,0,p1,kpx,0,kx,0)
	for y,kpy in enumerate(ky):


		Aint=[kpy]
		av = simulation(G,R,A,C,D,Aint,0,par,width,dx,time,dt,"TSXLT",time,False)
		#plot2D_simple(av,-1,axs,ki,1,n)
		plot2D_simple_tgh(av,-1,axs,x,1+y)
		#plot2D_kymograph_tgh(av,int(width/2), axs, x,1+y)
		niceaxis2(axs,x,y+1,p1,kpx,kpy,kx,ky)

plt.savefig("Figure/"+name+".pdf",dpi=300)
plt.show()

stop


#plot_crossection(av,-1,int(width/2),axs,5,1)
#plot_crosstime(av,int(width/2),axs,5,0)
'''




'''
====================================================================================================================================

TURING INSTABILTIY HERE

still in progress

======================================================================================================================================
'''

'''
filename="MAP_TuringInstability/type1"
p,df=loadTI(filename,parlist)

par_plot(df,parlist)
plt.savefig("Figures/"+"type1_"+"TI_Parplot.pdf",dpi=300)
plt.show()


filename="MAP_TuringInstability/type3"
p,df=load(filename,parlist)
par_plot(df,parlist)
plt.savefig("Figures/"+"type3_"+"TI_Parplot.pdf",dpi=300)
plt.show()

'''
'''
filename="MAP_TuringInstability/type1"
p,df=load(filename,parlist)
p0=pars_to_dict(p[8],parlist) #
#p0
par=addfixedpar(p0)


ttype, e2=getTuringinstability(par,200)
A0=np.zeros(1)
e2 = TuringInstability(A0,par,200)
print(e2.shape)
plt.axhline(y=0., xmin=0, xmax=200, c="k")
plt.plot(e2[0,0,0,:,:])
plt.yscale('symlog')
plt.show()

stop

'''





##########################################################################


'''
ttype, e2=getTuringinstability(par,200)
idx=np.argwhere(e2>0)[0]
ep=e2[idx[0],:,idx[2]]
axs[0,0].plot(ep,'r')
axs[0,0].set_yscale('symlog')
axs[0,0].axhline(y=0., xmin=0, xmax=200, c="blue")

k=np.linspace(1,5,20)
A0=np.logspace(-3,0,20)
A0[0]=10e-100
bifu_plot_tgh(par,A0,axs,0,1)

#A0[0]=0#np.logspace(-3,0,20)

#monostable_plot(A0,p0,axs,0,1,k,'K_GREEN2')

Tmatrx=[1,1,0]
Nmatrx=[0.,0.,0.,0.]
size=10 #mm 10-1
dx=0.1
par['D_ahl']=1.00
width=int(size/dx)
time=10
dt=0.01

A0=np.zeros(1)
ss=findss(A0,[0],par)
e = getEigen(ss,A0,par)

Vmatrx=[ss[0,0,0,0],ss[0,0,0,1],ss[0,0,0,2],1]
G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
av = simulation(G,R,A,C,D,[0],0,par,width,dx,time,dt,"TSXLT",time,False)
plot2D_simple_tgh(av,-1,axs,1,1)
plot2D_kymograph_tgh(av,int(width/2), axs, 2,1)

Vmatrx=[1,0.2,1,1]

Vmatrx=[ss[0,0,2,0],ss[0,0,2,1],ss[0,0,2,2],1]

G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
av = simulation(G,R,A,C,D,[0],0,par,width,dx,time,dt,"TSXLT",time,False)
plot2D_simple_tgh(av,-1,axs,1,0)	
plot2D_kymograph_tgh(av,int(width/2), axs, 2,0)


plt.show()




k=np.linspace(1,5,20)
A0=np.logspace(-3,0,20)
A0[0]=0#np.logspace(-3,0,20)

monostable_plot(A0,p0,axs,0,0,k,'K_GREEN')
plt.show()







stop


p=p[0:50]

A0=np.zeros(1)
e= TuringInstability(A0,p0,10)
#print(e)

ttype, e2=getTuringinstability(p0,10)
idx=np.argwhere(e2>0)[0]
ep=e2[idx[0],:,idx[2]]


s=int(np.round((np.sqrt(len(p))+.5)))
print(s)
fig, axs = plt.subplots(s,s,constrained_layout=True,figsize=(s,s))
x,y = 0, 0
for pi in p:
	p0=pars_to_dict(pi,parlist)
	p0=addfixedpar(p0)
	ttype, e2=getTuringinstability(p0,200)
	idx=np.argwhere(e2>0)[0]
	ep=e2[idx[0],:,idx[2]]
	axs[x,y].plot(ep,'r')
	axs[x,y].set_yscale('symlog')
	axs[x,y].axhline(y=0., xmin=0, xmax=200, c="blue")

	x+=1
	if x>s-1:
		x=0
		y+=1

plt.show()


#try to map max high

filename="MAP_TuringInstability/type1"
p,df=load(filename,parlist)
mm_tot=[]
for pi,pp in enumerate(p):
	p0=pars_to_dict(pp,parlist)
	p0=addfixedpar(p0)
	ttype, e2=getTuringinstability(p0,200)
	mm=np.nanmax(e2.real)
	mm_tot.append(mm)

par_plot(df,parlist,mm_tot)
plt.show()







stop

'''



#===============================================================================================================================================================


#Fit analysis 


#======================================================================================================================================================================



filename="FIT010_TSLT_gated"
modeltype="TSLT"
data="data_gated.txt"
datafile = 'data/'+modeltype + '/' +data


#n=['final','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18']
#n=['7']

#sys.path.insert(0, '/users/ibarbier/AC-DC/'+filename)
sys.path.insert(0, 'C:/Users/Administrator/Desktop/Modeling/TSReactionDiffusion/'+filename)
import model_fit as meq
parlist=meq.parlist

n=['43']
#n=['26']


gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
p, df= load(n[0],filename,meq.parlist)
p0=p[0]



p0['D_ahl']=0.1#0.010
p0['beta_ahl']=0#-1.5

p0['delta_ahl']=1
p0['K_ahl']=0
p0['n_ahl']=1.00


'''
#PARPLOT
par_plot(df,parlist)
plt.savefig(filename+"/plot/"+str(n[0])+"_parspace.pdf", dpi=300)
plt.show()
'''

'''
#COMPARE PLOT

#compare_plot([p[0]],filename,meq,datafile,modeltype,lw=1.5)
compare_plot(p,filename,meq,datafile,modeltype,lw=.5)
plt.savefig(filename+"/plot/"+str(n[0])+'_plot.pdf', bbox_inches='tight',dpi=300)
plt.show()
'''

'''
#BIFUPLOT
gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
fig, axs = plt.subplots(2,2,constrained_layout=True)
bifu_2Dplot_IPTG(AHL,p0,axs,0,0,IPTG,p,meq,"TSXLT")
size=50
AHL=np.logspace(-5,1,size)
IPTG=np.logspace(-2,1,size)
bifu_2Dplot_IPTG(AHL,p0,axs,1,0,IPTG,p,meq,"TSXLT")
plt.savefig(filename+"/plot/"+str(n[0])+'_bifuDiag.pdf', bbox_inches='tight',dpi=300)
plt.show()
'''

#TURING Instability







'''
gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
nPar=9
k1=np.linspace(-3,0,nPar)
k2=np.linspace(-2,2,nPar)
par1="beta_ahl"
par2="K_ahl"
fig, axs = plt.subplots(nPar,nPar)#constrained_layout=True)
for x,px in enumerate(k1):
	for y,py in enumerate(k2):
		p0[par1]=px
		p0[par2]=py
		bifu_2Dplot_IPTG(AHL,p0,axs,x,y,IPTG,p,meq,"TSXLT")
		niceaxis(axs,x,y,px,py,k1,k2,8)
plt.show()
'''


#colonies simulation

########################################



A0=np.logspace(-6,0,100)
I=np.ones(1)*0.25

ss=meq.findss(A0,I,p0,"TSXLT")
gstate=np.copy(ss[0,0,1])
gstate[0]=gstate[0]-10**p0['basal_green']
#gstate[2]=gstate[2]-10**p0['basal_red']

#rstate=np.copy(ss[-1,-1,0])
#rstate[0]=rstate[0]-10**p0['basal_green']
#rstate[1]=rstate[1]-10**p0['basal_red']


fig, axs = plt.subplots(2,3)#constrained_layout=True)


bifu_plot_fit(p0,A0,I,axs,0,0,meq,"TSXLT")
plt.show()
A0=np.zeros(1)


Tmatrx=[1,1,0] #which ?
Nmatrx=[0.,0.,0.,0.] #noise
Vmatrx=[800,0,0,1] #initial values
size= 16
dx=0.05
time=4
dt=0.001
width=int(size/dx)

fig, axs = plt.subplots(2,2)#constrained_layout=True)
G,R,A,C,D = init_colony(Vmatrx,width)
av = simulation(G,R,A,C,D,A0,I,p0,width,dx,time,dt,"TSXLT_fit",time,True)
plot2D_simple_tgh(av,-1,axs,0,0)
plot2D_kymograph_tgh(av,int(width/2), axs, 0, 1)
#plot_crossection(av,-1,int(width/2), axs, 0, 1)
plot_crossection_diffusion(av,-1,int(width/2), axs, 1, 1)
plot_crosstime(av,int(width/2), axs, 1, 0)


plt.show()


stop

#2D simulation
Tmatrx=[1,1,0]
Nmatrx=[0.,0.,0.,0.]


width=int(size/dx)
time=200
dt=0.1

A0=AHL
I=IPTG[:,np.newaxis]

gstate=np.copy(ss[0,0,0])
gstate[0]=gstate[0]-10**p0['basal_green']
gstate[1]=gstate[1]-10**p0['basal_red']

rstate=np.copy(ss[-1,-1,0])
rstate[0]=rstate[0]-10**p0['basal_green']
rstate[1]=rstate[1]-10**p0['basal_red']

fig, axs = plt.subplots(2,2,constrained_layout=True)

Vmatrx=[gstate[0],gstate[1],0,1]
G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
av = simulation(G,R,A,C,D,A0,I,p0,width,dx,time,dt,"TSXLT_fit",time,False)
plot2D_simple_tgh(av,-1,axs,0,0)

Vmatrx=[rstate[0],rstate[1],0,1]
G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
av = simulation(G,R,A,C,D,A0,I,p0,width,dx,time,dt,"TSXLT_fit",time,False)
plot2D_simple_tgh(av,-1,axs,1,0)
plt.show()

stop







#===============================================================================================================================================================


#delete?

#======================================================================================================================================================================


##################################################################333

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
				C=np.ones((width,width))
				D=np.ones((width,width))

			elif v==1:
				G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
				R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				C=np.ones((width,width))
				D=np.ones((width,width))
			elif v==2:
				G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
				C=np.ones((width,width))
				D=np.ones((width,width))
			elif v==3:
				G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
				A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
			elif v==4:
				G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
				R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				C=np.ones((width,width))
				D=np.ones((width,width))
			elif v==5:
				G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
				C=np.ones((width,width))
				D=np.ones((width,width))

			av=simulation(G,R,A,C,D,A0,0,par,width,dx,time,dt,model,tt=time)

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
					C=np.ones((width,width))
					D=np.ones((width,width))					

				elif v==1:
					G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))*0
					R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					C=np.ones((width,width))
					D=np.ones((width,width))
				elif v==2:
					G=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					R=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					A=np.random.normal(loc=0.8, scale=0.2, size=(width,width))
					C=np.ones((width,width))
					D=np.ones((width,width))
				elif v==3:
					G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
					A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					C=np.ones((width,width))
					D=np.ones((width,width))
				elif v==4:
					G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))*0
					R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					C=np.ones((width,width))
					D=np.ones((width,width))
				elif v==5:
					G=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					R=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					A=np.random.normal(loc=0.1, scale=0.1, size=(width,width))
					C=np.ones((width,width))
					D=np.ones((width,width))

				av=simulation(G,R,A,C,D,A0,0,par,width,dx,time,dt,model,tt=time)

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






par['K_ahl_green']=0.#k2[15]#[5]
par['K_ahl_red']= 0.

par['K_GREEN']=1#k[8]#[3]
par['K_RED']= 1 #ii

par['delta_red']=1.
par['delta_green']=1
par['D_ahl']=0.1
'''

'''
A0=np.logspace(-3,0,width)
bifu_plot(A0,par)
plot2D(par)
'''
#plot2D(par)


#A0=np.logspace(-2,0,10)
#k = np.linspace(0.,2,size) #for K_RED/GREEN
#screen_2par(par,'K_GREEN','K_RED',k,k)
k = np.linspace(-1,1,nPar) #for K_RED/GREEN



plot_all(par,model="TSXLT",k=k,parm="K_ahl_green", i='w')

par=par_turing
#plot1D(par,model="Turing")

plot2D(par,model="Turing")


'''
I=0
A0=np.logspace(-3,0,width)
k = np.linspace(0,2,size) #for K_RED/GREEN
k2 = np.linspace(-1,1,size) #for K_ahl
k3 = np.linspace(1,3,size) #for K_ahl

bifu_plot(A0,par)
plot_eigen(par)
plot1D(par,endpoint=True)

A0=np.ones(1)*10**-0.5

ss=findss(A0,[0],par)

G=np.ones(1)*ss[-1,:,2,0]
R=np.ones(1)*ss[-1,:,2,1]
A=np.ones(1)*ss[-1,:,2,2]

av=simulation(G-10e-1,R,A,A0,0,par,dt,model="TSXLT",tt=time)

plt.plot(av[:,0])
#plt.plot(av[:,1])
#plt.plot(av[:,2])
plt.show()


'''

'''
print("-----------------------------------------------")


A0=np.ones(1)*10**0

ss=findss(A0,[0],par)
print(ss)

G=np.ones(1)*ss[-1,:,0,0]
R=np.ones(1)*ss[-1,:,0,1]
A=np.ones(1)*ss[-1,:,0,2]

print("SS: ")
print(G,R,A)

sse=getEigen(ss,A0,par)
print("eigen: ")
print( sse)

g,r,a = model_TSXLT(G,R,A,A0,I,par)


print("dx: ")
print(g,r,a)

av=simulation(G,R,A,A0,0,par,dt,model="TSXLT",tt=time)
print("sim: " )
print(av[-1,0],av[-1,1],av[-1,2])


ss2=np.array([[[(av[-1,0][0],av[-1,1][0],av[-1,2][0])]]])
sse2=getEigen(ss2,A0,par)

print("eigen: ")

print(sse2)

print("------------------------------")

print([av[-1,0][0]==G,av[-1,1][0]==R,av[-1,2][0]==A])

plt.plot(av[:,0])
#plt.plot(av[:,1])
#plt.plot(av[:,2])
plt.show()
print("------------------------------")

'''








