
'''
Here are the script to run stuff

written by Içvara
'''
from model import *
from diffusion import *
from torun import *

from grow import *
from  easyplot import *
from TuringInstability import*
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import sys






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

    'D_ahl':0,#0.1,

    'K_IPTG':1
    }

'''
size=2 #mm 10-1
dx=0.05
width=int(size/dx)

#dt= 0.9 * (dx ** 2) / (2 * D)

#ar0['D_ahl']=.01
dx_diffusion=0.05
time=16
dt=0.005
'''



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

size=4 #cm 10-1
#dt=0.005
dt=0.1
size=4
dx=0.05
width=int(size/dx)
pos=int(width/2)
unitX= np.arange(0,width)*dx
time=20
par_growth['H_growth'] = 0.05 #0.05# #  cm/h ?
par_growth['D_growth'] =  0.5 #  cm/h ?

######## 1. TEST DIFFUSION ###########


#test_diffusion("1")
#stop



######### TEST GROWTH  stat + heredity####################

'''
par=par0
par=par_rep

par_growth['H_growth'] = 0.05 #0.05# #  cm/h ?
par_growth['D_growth'] =  0.5 #  cm/h ?
A0=np.zeros(1)
#par_growth['max_density']=10


Tmatrx=[1,1,1] # G,R,AHL  tell which one to randomize
Vmatrx=[0,1.,0,1] # same order, the init values
Nmatrx=[0.,0.,0.,0.] #the noise itensity G,R,A,D

#G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
G,R,A,C,D = init_colony(Vmatrx,width)

'''
'''
name ="split/"
G[pos,pos]=0
R[pos,pos]=1
C[pos,pos]=1
D[pos,pos]=1

G[pos-1,pos-1]=0
R[pos-1,pos-1]=1
C[pos-1,pos-1]=1
D[pos-1,pos-1]=1

G[pos,pos-1]=1
R[pos,pos-1]=0
C[pos,pos-1]=1
D[pos,pos-1]=1

G[pos-1,pos]=1
R[pos-1,pos]=0
C[pos-1,pos]=1
D[pos-1,pos]=1
'''
'''

#name='TSXL_R_lawn2'
fig, axs = plt.subplots(2,2,constrained_layout=True)

#av = simulation(G,R,A,C,D,A0,0,par,width,dx,time,dt,"TSXLT",time,False)
av = simulation(G,R,A,C,D,A0,0,par,width,dx,time,dt,"Repressilator",time,True)

#plot2D_simple_tgh(av,-1,axs,0,0)
plot2D_simple_tgh3(av,-1,axs,0,0)
axs[0,1].imshow(av[-1,5,:,:],cmap="viridis") #, norm=LogNorm())
axs[1,0].plot(av[-1,0,pos,pos].T,color="g")
axs[1,0].plot(av[-1,1,pos,pos].T,color="r")
axs[1,0].plot(av[-1,3,pos,pos].T,color="b")


axs[1,1].plot(unitX,av[:,5,:,pos].T,color="m")
axs[1,1].set_yscale("log")
plt.show()
plt.close()


name ="rep_sta/"
for t in np.arange(av.shape[0]):
	fig, axs = plt.subplots(2,2,constrained_layout=True)
	#plot2D_simple_tgh(av,t,axs,0,0)
	plot2D_simple_tgh3(av,t,axs,0,0)

	axs[0,1].imshow(av[t,5,:,:],cmap="viridis") #, norm=LogNorm())
	#axs[1,0].plot(unitX,av[t,4,:,pos].T,color='k')
	axs[1,0].plot(unitX,av[t,0,:,pos].T,color='g')
	axs[1,0].plot(unitX,av[t,1,:,pos].T,color='r')
	axs[1,0].plot(unitX,av[t,3,:,pos].T,color='b')

	axs[1,1].plot(unitX,av[t,5,:,pos].T,color="m")
	axs[1,1].set_yscale("log")
	axs[1,1].set_ylim(1,100)
	plt.savefig("Figures/timelapse/"+name+ str(t) +".png",dpi=300)
	plt.close()

stop
'''





###################


'''
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

		#plot2D_simple_tgh(av,-1,axs,x,1+y)
		#plot2D_kymograph_tgh(av,int(width/2), axs, x,1+y)
#plot_crossection(av,-1,int(width/2),axs,5,1)
#plot_crosstime(av,int(width/2),axs,5,0)
'''





'''
====================================================================================================================================

TURING INSTABILTIY HERE

still in progress

======================================================================================================================================
'''

############3
# par space
###################

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

############3
# Eigens map
###################

'''
filename="MAP_TuringInstability/type1"
p,df=loadTI(filename,parlist)
#p0=pars_to_dict(p[8],parlist) #
p0=p[8]
print(p0)
par=addfixedpar(p0)
ttype, e2=getTuringinstability(par,200)
A0=np.zeros(1)
e2 = TuringInstability(A0,par,200)
plt.axhline(y=0., xmin=0, xmax=200, c="k")
plt.plot(e2[0,0,0,:,:])
plt.yscale('symlog')
plt.savefig("Figures/"+"type1_"+"plot.pdf",dpi=300)
plt.show()
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

n=['26']
n=['48']


gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
p, df= load(n[0],filename,meq.parlist)
pMAP=pars_to_dict(getMAP(df.to_numpy().T),parlist)
p0=p[0]
pmedian=pars_to_dict(np.array(df.median().tolist()[:-1]),parlist)

print(pMAP)
print(p[239])


stop

#p0['delta_red']=1
#p0['delta_green']=1


'''
pMAP['D_ahl']=1#0.010
pMAP['beta_ahl']=-100.5#-1.5
pMAP['delta_ahl']=1
pMAP['K_ahl']=0
pMAP['n_ahl']=1.00
'''


###########
#PARPLOT
###########
'''
par_plot(df,parlist)
plt.savefig(filename+"/plot/"+str(n[0])+"_parspace.pdf", dpi=300)
plt.show()
'''

###########
#COMPARE PLOT
###########

'''
compare_plot([p[0]],filename,meq,datafile,modeltype,lw=1.5)
plt.savefig(filename+"/plot/"+str(n[0])+'p0_plot.pdf', bbox_inches='tight',dpi=300)
plt.show()
compare_plot([pMAP],filename,meq,datafile,modeltype,lw=1.5)
plt.savefig(filename+"/plot/"+str(n[0])+'pMAP_plot.pdf', bbox_inches='tight',dpi=300)
plt.show()

compare_plot(p,filename,meq,datafile,modeltype,lw=.5)
plt.savefig(filename+"/plot/"+str(n[0])+'_plot.pdf', bbox_inches='tight',dpi=300)
plt.show()
'''


###########
#COMPARE PLOT from TSLT to TSXLT
###########
#doesnt' work yet
#modeltype="TSXLT"
#datafile = 'data/'+modeltype + '/' +data

#p0['K_RED']=-0.5sss
#compare_plot([pMAP],filename,meq,datafile,modeltype,lw=1.5)
#plt.savefig(filename+"/plot/"+str(n[0])+'_pTEST_plot.pdf', bbox_inches='tight',dpi=300)
#plt.show()
###########
#BIFUPLOT 2D
###########

'''
gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
fig, axs = plt.subplots(2,2,constrained_layout=True)
bifu_2Dplot_IPTG(AHL,pMAP,axs,0,0,IPTG,p,meq,"TSXLT",1)
size=200
AHL=np.logspace(-6,0,size)
AHL=np.logspace(-4,0,size)


IPTG=np.logspace(-2.5,1,size)
IPTG=np.logspace(-5,1,size)

bifu_2Dplot_IPTG(AHL,pMAP,axs,1,0,IPTG,p,meq,"TSXLT",1)

AHL=np.logspace(-6,0,size)
AHL=np.logspace(-4,0,size)

IPTG=np.logspace(-5,0,size)
bifu_2Dplot_IPTG(AHL,pMAP,axs,1,1,IPTG,p,meq,"TSXLT",1)

plt.savefig(filename+"/plot/"+'pMAP_bifuDiag.pdf', bbox_inches='tight',dpi=300)
plt.show()

stop
'''
###########
#BIFUPLOT 1D
###########

'''
gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
size=6
fig, axs = plt.subplots(size,2,constrained_layout=True, figsize=(4,size))
AHL=np.logspace(-6,0,200)
#IPTG=np.logspace(-2.5,0,size)
#size=20
#IPTG=np.logspace(-0.2,0.2,size)
bifu_plot_fit(pMAP,AHL,IPTG,axs,0,0,meq,"TSXLT")
plt.savefig(filename+"/plot/"+str(n[0])+'_bifucurve2.pdf', bbox_inches='tight',dpi=300)
plt.show()
'''

###########
#TURING Instability
###########

'''
gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
TI_plot_fit(AHL,IPTG,pMAP,meq)
plt.savefig(filename+"/plot/"+str(n[0])+'_pMAPTI.pdf', bbox_inches='tight',dpi=300)
plt.show()
'''



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


n=['26']
#n=['48']
p, df= load(n[0],filename,meq.parlist)
pMAP=pars_to_dict(getMAP(df.to_numpy().T),parlist)
pp=p[239]
#pp=pMAP


#pMAP['D_ahl']=1#0.010
'''
pMAP['beta_ahl']=-100.5#-1.5
pMAP['delta_ahl']=1
pMAP['K_ahl']=0
pMAP['n_ahl']=1.00
'''



#dt=0.005
dt=0.005
size=4
dx=0.05
width=int(size/dx)
pos=int(width/2)
unitX= np.arange(0,width)*dx
time=20#20#20
par_growth['H_growth'] = 0.05 #0.05# #  cm/h ?
par_growth['D_growth'] =  0.5 #  cm/h ?
par_growth['D_growth'] =  1 #  cm/h ?



#dt=0.05
'''
par=pMAP.copy()
par['delta_green']=1
par['delta_red']=1
par['delta_ahl']=1

par['beta_red']=1
par['alpha_red']=3
par['K_ahl_green']=1
par['K_ahl_red']=1
par['beta_green']=1
par['alpha_green']=3

par['basal_green']=0.
par['basal_red']=0.


#par=par_TSLT_mushroom
par['K_GREEN']=2 # here == IPTG
par['K_RED']=1

par['n_GREEN']=2
par['n_RED']=2
'''

#par['D_ahl']=size/10/time #mm/h
IPTG=np.logspace(-2.5,1,6) #6


gg,gr,rr,rg=meq.Get_data(datafile)
AHL=gg.index.values
IPTG=gg.columns.values
A0= np.ones(1)*.0
A0_bifu=np.logspace(-5,1,100)

AHL=np.logspace(-6,0,6)
#AHL[0]=0
IPTG=np.logspace(-2.5,1,6)
I0=np.ones(1)*0.5
ss=meq.findss(AHL,IPTG,pp,modeltype)
print(ss)

Tmatrx=[1,1,0]
state=0
#Vmatrx=[ss[-1,0,state,0]-10**pp['basal_green'],ss[-1,0,state,2],0,1]

Vmatrx=[ss[0,0,state,0]-10**pp['basal_green'],ss[0,0,state,2],ss[0,0,state,3],1]

Vmatrx2=[ss[-1,-1,state,0]-10**pp['basal_green'],ss[-1,-1,state,2],ss[0,0,state,3],1]

#Vmatrx=[ss[0,0,state,0]-10**pp['basal_green'],ss[0,0,state,2],ss[0,0,state,3],1]

print(Vmatrx)
print(Vmatrx2)

Vmatrx=[ss[0,0,state,0]-10**pp['basal_green'],ss[0,0,state,2],ss[0,0,state,3],1]



Nmatrx=[0,0,0,0]

G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)

G[40,40]=Vmatrx2[0]
R[40,40]=Vmatrx2[1]
A[40,40]=Vmatrx2[2]
#G,R,A,C,D = init_colony(Vmatrx,width)





name='lawn_TSXLT_dif'
pp['D_ahl']=0.#0.1 #cm/h 2.5 10**-3 cm2/h (Basu et al, 2005; Miyamoto et al,2018). 



fig, axs = plt.subplots(6,6,constrained_layout=True)

print(pp)
#pp['beta_ahl']=0

#bifu_2Dplot_IPTG(AHL,pp,axs,0,0,IPTG,p,meq,"TSXLT",1)
IPTG_bifu=np.logspace(-2.5,1,100)

IPTG=np.logspace(-2.5,1,6) #6
#IPTG= 
#IPTG=gg.columns.values
#IPTG[0]=10**-2.5
#IPTG=IPTG[-2]
#for i,I in enumerate([IPTG[0]]):

for i,I in enumerate(IPTG):
	bifu_plot_fit_tgh2(pp,A0,IPTG_bifu,axs,i,0,meq,"TSXLT")





#bifu_plot_fit_tgh(pp,A0_bifu,I0,axs,1,0,meq,"TSXLT")

	axs[i,0].axvline(x=np.log10(I))

#col sim
	for j,n in enumerate([0]):#[0,0.5,1,2,4]):

		av = simulation(G,R,A,C,D,A0,I,pp,width,dx,time,dt,"TSXLT_fit",time,n,False,True,meq=meq)
		#av2 = simulation(G,R,A,C,D,A0,I,pp,width,dx,time,dt,"TSXLT_fit",time,n,True,meq=meq)
		ngraph=1
		plot2D_simple_tgh(av,-1,axs,i,(j*ngraph+1))
		plot2D_kymograph_tgh(av,int(width/2),axs,i,(j*ngraph+2))

		axs[i,(j*ngraph+3)].plot(av[-1,0,int(width/2),:],'g')
		axs[i,(j*ngraph+3)].plot(av[-1,3,int(width/2),:],'b')
		axs[i,(j*ngraph+3)].plot(av[-1,1,int(width/2),:],'r')

		'''
		plot2D_simple_tgh(av2,-1,axs,i,(j+j*ngraph+3))
		axs[i,(j+j*ngraph+4)].plot(av2[-1,0,int(width/2),:],'g')
		axs[i,(j+j*ngraph+4)].plot(av2[-1,3,int(width/2),:],'b')
		axs[i,(j+j*ngraph+4)].plot(av2[-1,1,int(width/2),:],'r')
		'''



plt.savefig("Figures/"+name+".pdf",dpi=300)



plt.show()


stop


#############################33







A0=np.logspace(-6,0,100)
I=np.ones(1)*0.25

ss=meq.findss(A0,I,pMAP,"TSXLT")
gstate=np.copy(ss[0,0,1])
gstate[0]=gstate[0]-10**pMAP['basal_green']
#gstate[2]=gstate[2]-10**p0['basal_red']

#rstate=np.copy(ss[-1,-1,0])
#rstate[0]=rstate[0]-10**p0['basal_green']
#rstate[1]=rstate[1]-10**p0['basal_red']


fig, axs = plt.subplots(2,3)#constrained_layout=True)


bifu_plot_fit(pMAP,A0,I,axs,0,0,meq,"TSXLT")
plt.show()
A0=np.zeros(1)


Tmatrx=[1,1,0] #which ?
Nmatrx=[0.,0.,0.,0.] #noise
Vmatrx=[800,0,0,1] #initial values
size= 16
dx=0.05
time=1
dt=0.001
width=int(size/dx)

fig, axs = plt.subplots(2,2)#constrained_layout=True)
G,R,A,C,D = init_colony(Vmatrx,width)
av = simulation(G,R,A,C,D,A0,I,pMAP,width,dx,time,dt,"TSXLT_fit",time,True,meq=meq)
plot2D_simple_tgh(av,-1,axs,0,0)
plot2D_kymograph_tgh(av,int(width/2), axs, 0, 1)
#plot_crossection(av,-1,int(width/2), axs, 0, 1)
plot_crossection_diffusion(av,-1,int(width/2), axs, 1, 1)
plot_crosstime(av,int(width/2), axs, 1, 0)


plt.show()



#2D simulation
Tmatrx=[1,1,0]
Nmatrx=[0.,0.,0.,0.]


width=int(size/dx)
time=200
dt=0.1

A0=AHL
I=IPTG[:,np.newaxis]

gstate=np.copy(ss[0,0,0])
gstate[0]=gstate[0]-10**pMAP['basal_green']
gstate[1]=gstate[1]-10**pMAP['basal_red']

rstate=np.copy(ss[-1,-1,0])
rstate[0]=rstate[0]-10**pMAP['basal_green']
rstate[1]=rstate[1]-10**pMAP['basal_red']

fig, axs = plt.subplots(2,2,constrained_layout=True)

Vmatrx=[gstate[0],gstate[1],0,1]
G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
av = simulation(G,R,A,C,D,A0,I,pMAP,width,dx,time,dt,"TSXLT_fit",time,False)
plot2D_simple_tgh(av,-1,axs,0,0)

Vmatrx=[rstate[0],rstate[1],0,1]
G,R,A,C,D = init_grid(Vmatrx,Nmatrx,width,Tmatrx)
av = simulation(G,R,A,C,D,A0,I,pMAP,width,dx,time,dt,"TSXLT_fit",time,False)
plot2D_simple_tgh(av,-1,axs,1,0)
plt.show()

stop





