'''
Here are the different equation for TS system and Turing system

There is the equation for simulation and the function to find steady state

written by Içvara
'''

import numpy as np
from scipy.signal import argrelextrema
from scipy import optimize
from scipy.optimize import brentq
import matplotlib.pyplot as plt


par_turing={

    'beta_green':1,
    'K_ahl_green':1.5,
    'n_ahl_green':2,
    'delta_green':1,


    'beta_red':1.00,
    'K_ahl_red':1,
    'n_ahl_red':2,
    'delta_red':1.,

    'K_RED':2.5,
    'n_RED':2,


    'D_red':1,
    'D_green':0.01

    }


par_rep={

    'beta_green':2,
    'delta_green':1,


    'beta_red':2.00,
    'delta_red':1.,

    'beta_blue':2.00,
    'delta_blue':1.,

    'K_RED':2,
    'n_RED':2,

    'K_GREEN':2,
    'n_GREEN':2,

    'K_BLUE':2,
    'n_BLUE':2,



    }

#sometimes with power it replace it by 0 , but make the script buggy
minabs=10e-300


def model_turing(GREENi,REDi,par):

    GREEN = 0.01 + par['beta_green']*np.power(GREENi*10**par['K_ahl_green'],par['n_ahl_green'])/(1+np.power(GREENi*10**par['K_ahl_green'],par['n_ahl_green']))
    GREEN = GREEN / (1 + np.power(REDi*10**par['K_RED'],par['n_RED']))
    GREEN = GREEN - par['delta_green']*GREENi  
    

    RED = (par['beta_red']*np.power(GREENi*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(GREENi*10**par['K_ahl_red'],par['n_ahl_red']))
   # RED = RED / (1 + np.power(GREENi*10**par['K_GREEN'],par['n_GREEN']))
    RED = RED - par['delta_red']*REDi 

    return GREEN,RED


def model_TSXLT(GREENi,REDi,AHLi,A0,IPTG,par):
    '''
    maxmax=1e50
    minmin=-1e50
    GREENi[GREENi>maxmax]=maxmax
    REDi[REDi>maxmax]=maxmax
    AHLi[AHLi>maxmax]=maxmax

    GREENi[GREENi<minmin]=minmin
    REDi[REDi<minmin]=minmin
    AHLi[AHLi<minmin]=minmin
    '''

    AHLii =  AHLi + A0

    GREEN = par['alpha_green'] + par['beta_green']*np.power(AHLii*10**par['K_ahl_green'],par['n_ahl_green'])/(1+np.power(AHLii*10**par['K_ahl_green'],par['n_ahl_green']))
    #GREEN[np.isnan(GREEN)] = minabs
    GREEN = GREEN / (1 + np.power(REDi*10**par['K_RED'],par['n_RED']))
    GREEN = GREEN - par['delta_green']*GREENi  
    free_GREENi= GREENi / ( 1+ 10**par['K_IPTG']*IPTG)

    RED = par['alpha_red']+(par['beta_red']*np.power(AHLii*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(AHLii*10**par['K_ahl_red'],par['n_ahl_red']))
    #RED[np.isnan(RED)] = minabs
    RED = RED / (1 + np.power(free_GREENi*10**par['K_GREEN'],par['n_GREEN']))
    RED = RED - par['delta_red']*REDi 

    AHL = (par['beta_ahl']*np.power(GREENi*10**par['K_ahl'],par['n_ahl']))/(1+np.power(GREENi*10**par['K_ahl'],par['n_ahl']))
    AHL = AHL - par['delta_ahl']*AHLi 

    return GREEN,RED,AHL#,RED_Giac



def model_Repressilator(GREENi, REDi, BLUEi,par):

    GREEN = par['beta_green'] / (1 + np.power(BLUEi*10**par['K_BLUE'],par['n_BLUE']))
    GREEN = GREEN - par['delta_green']*GREENi 

    RED = par['beta_red'] / (1 + np.power(GREENi*10**par['K_GREEN'],par['n_GREEN']))
    RED = RED - par['delta_red']*REDi 

    BLUE = par['beta_blue'] / (1 + np.power(REDi*10**par['K_RED'],par['n_RED']))
    BLUE = BLUE - par['delta_blue']*BLUEi 

    return GREEN,RED,BLUE






   
def solvedfunction(Gi,Ai,I,par,model="TSXLT"):
    #rewrite the system equation to have only one unknow and to be call with scipy.optimze.brentq
    #the output give a function where when the line reach 0 are a steady states
    if model != "TSXLT":
        return 0

    Gf = Gi / ( 1+ 10**par['K_IPTG']*I)

    A = (par['beta_ahl']*np.power(Gi*10**par['K_ahl'],par['n_ahl']))/(1+np.power(Gi*10**par['K_ahl'],par['n_ahl']))
    A = A / par['delta_ahl']

    #A= AHL + Diffusion(Ai)
    Aii = A  + Ai


    R = par['alpha_red'] + par['beta_red']*np.power(Aii*10**par['K_ahl_red'],par['n_ahl_red'])/(1+np.power(Aii*10**par['K_ahl_red'],par['n_ahl_red']))
    R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN']))  
    R = ( R ) / par['delta_red']  #################


    G = par['alpha_green'] +( par['beta_green']*np.power(Aii*10**par['K_ahl_green'],par['n_ahl_green']))/(1+np.power(Aii*10**par['K_ahl_green'],par['n_ahl_green'])) 
    G = G / (1 + np.power(R*10**par['K_RED'],par['n_RED'])) #+ 10**par['basal_green']
    G = (G ) / par['delta_green']

    func = G - Gi

    return func 


def findss(Ai,I,par,model="TSXLT"):
    if model != "TSXLT":
        return 0
    ss=[]
    nNode=3 # number of nodes : X,Y,Z
    nStstate= 5
    nAHL= len(Ai)
    nIPTG=len(I)
    ss=np.ones((nAHL,nIPTG,nStstate,nNode))*np.nan  
    for ai,a in enumerate(Ai):
        for iptgi,iptg in enumerate(I):
            Gi=np.logspace(-10,5,1000,base=10)
            f=solvedfunction(Gi,a,iptg,par)
            x=f[1:-1]*f[0:-2] #when the output give <0, where is a change in sign, meaning 0 is crossed
            index=np.where(x<0)
            for it,i in enumerate(index[0]):
                G = brentq(solvedfunction, Gi[i], Gi[i+1],args=(a,iptg,par)) #find the value of AHL at 0
                Gf = G / ( 1+ 10**par['K_IPTG']*iptg)

                A = (par['beta_ahl']*np.power(G*10**par['K_ahl'],par['n_ahl']))/(1+np.power(G*10**par['K_ahl'],par['n_ahl']))
                A = A / par['delta_ahl']
                Aii= A + a
                R =par['alpha_red'] + par['beta_red']*np.power(Aii*10**par['K_ahl_red'],par['n_ahl_red'])/(1+np.power(Aii*10**par['K_ahl_red'],par['n_ahl_red']))
                R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN']))  
                R = ( R ) / par['delta_red']  #################

                ss[ai,iptgi,it]=np.array([G,R,A])
             #   ss[ai,iptgi,it]=np.array([G+par['basal_green'],R+par['basal_red']])
    return ss


def jacobianMatrix(ss,par,A0,model="TSXLT"):
    #need ss in [AHL,fluo] dimension

    JM=np.ones((ss.shape[0],ss.shape[1],ss.shape[2],3,3))*np.nan 
    I=0
    if model=="TSXLT":
        G=np.copy(ss[:,I,:,0])
        R=np.copy(ss[:,I,:,1])
        A=np.copy(ss[:,I,:,2])

        A=A+A0[:,np.newaxis]

        dGdg = - par['delta_green'] *G/G

        dGdr = -(((par['beta_green']*np.power((10**par['K_ahl_green']*A),par['n_ahl_green']))/(1+np.power((10**par['K_ahl_green']*A),par['n_ahl_green']))+par['alpha_green'])*par['n_RED']*np.power((10**par['K_RED']*R),par['n_RED']))
        dGdr =  dGdr/(R*np.power((np.power((10**par['K_RED']*R),par['n_RED'])+1),2))

        dGda = par['beta_green']*par['n_ahl_green']*np.power((A*10**par['K_ahl_green']),par['n_ahl_green'])
        dGda= dGda/ ((np.power(R*10**par['K_RED'],par['n_RED'])+1)*A*np.power((np.power(A*10**par['K_ahl_green'],par['n_ahl_green'])+1),2))
           
        dRdr = - par['delta_red'] *R/R

        dRdg = -(((par['beta_red']*np.power((10**par['K_ahl_red']*A),par['n_ahl_red']))/(1+np.power((10**par['K_ahl_red']*A),par['n_ahl_red']))+par['alpha_red'])*par['n_GREEN']*np.power((10**par['K_GREEN']*G),par['n_GREEN']))
        dRdg =  dRdg/(G*np.power((np.power((10**par['K_GREEN']*G),par['n_GREEN'])+1),2))

        dRda = par['beta_red']*par['n_ahl_red']*np.power((A*10**par['K_ahl_red']),par['n_ahl_red'])
        dRda= dRda/ ((np.power(G*10**par['K_GREEN'],par['n_GREEN'])+1)*A*np.power((np.power(A*10**par['K_ahl_red'],par['n_ahl_red'])+1),2))

        dAdg = par['beta_ahl']*par['n_ahl']*np.power((G*10**par['K_ahl']),par['n_ahl'])
        dAdg = dAdg/ (G*np.power((np.power((10**par['K_ahl']*G),par['n_ahl'])+1),2))

        dAdr = 0 * A/A

        dAda =  - par['delta_ahl'] *A/A
           
        JM[:,0,:,0,0]=dGdg
        JM[:,0,:,0,1]=dGdr
        JM[:,0,:,0,2]=dGda

        JM[:,0,:,1,0]=dRdg
        JM[:,0,:,1,1]=dRdr
        JM[:,0,:,1,2]=dRda

        JM[:,0,:,2,0]=dAdg
        JM[:,0,:,2,1]=dAdr
        JM[:,0,:,2,2]=dAda
    else:
        JM=0   
    return JM


def approximateJacob(G,R,A,A0,I,par): 
    #allow to check if jacobian matrix derivate are correctly written
    delta=10e-5
    g,r,a = model_TSXLT(G,R,A,A0,I,par)
    #print(g,r,a)

    dgdg = (model_TSXLT(G+delta,R,A,A0,I,par)[0]-g)/delta
    dgdr = (model_TSXLT(G,R+delta,A,A0,I,par)[0]-g)/delta
    dgda = (model_TSXLT(G,R,A+delta,A0,I,par)[0]-g)/delta

    drdg = (model_TSXLT(G+delta,R,A,A0,I,par)[1]-r)/delta
    drdr = (model_TSXLT(G,R+delta,A,A0,I,par)[1]-r)/delta
    drda = (model_TSXLT(G,R,A+delta,A0,I,par)[1]-r)/delta

    dadg = (model_TSXLT(G+delta,R,A,A0,I,par)[2]-a)/delta
    dadr = (model_TSXLT(G,R+delta,A,A0,I,par)[2]-a)/delta
    dada = (model_TSXLT(G,R,A+delta,A0,I,par)[2]-a)/delta

    JM=np.ones((G.shape[0],3,3))*np.nan 

    JM[:,0,0]=dgdg
    JM[:,0,1]=dgdr
    JM[:,0,2]=dgda

    JM[:,1,0]=drdg
    JM[:,1,1]=drdr
    JM[:,1,2]=drda

    JM[:,2,0]=dadg
    JM[:,2,1]=dadr
    JM[:,2,2]=dada

    return JM

def getEigen(ss,A0,par,model="TSXLT"):
    if model !="TSXLT":
        return 0
    J=jacobianMatrix(ss,par,A0,model)
    J2=np.nan_to_num(J)
    eigvals, eigvecs =np.linalg.eig(J2)
    sse=eigvals
    return sse #, np.trace(A), np.linalg.det(A)


def TuringInstability(A0,par,n=100,model="TSXLT"):
    if model !="TSXLT":
        return 0
    I=0
    q=np.linspace(0,20,n)#100
    ss=findss(A0,[I],par)
    J=jacobianMatrix(ss,par,A0,model)
    JE=np.ones((len(A0),len([I]),5,n,3,3))*np.nan
    JE[:,:,:,:,:,:]=np.copy(J[:,:,:,np.newaxis,:,:])
    #JE=np.ones((len(A0),len([I]),5,3,3,20))*np.copy(J)

    eigens=np.ones((len(A0),len([I]),5,len(q),3))*np.nan


    if model =="TSXLT":
        JE[:,:,:,:,2,2]=JE[:,:,:,:,2,2]- q**2*par['D_ahl']
        JE2=np.nan_to_num(JE)
        eigvals, eigvecs =np.linalg.eig(JE2)
        sse=eigvals.real
        sse=np.sort(sse,4)
    return sse


'''
def bifu_plot(par,ID,inter):
    IPTG=np.logspace(-2,1,50)
    AHLi=np.logspace(-2,0,100)
    ss = findss(AHLi,IPTG,par)
    fig, axs = plt.subplots(3, 3)
    for mi,m in enumerate(['go','yo','go','yo','ro']):
                    axs[0,0].plot(ss[:,0,mi,0],m,markersize=1.)
                    axs[1,0].plot(ss[:,inter,mi,0],m,markersize=1.)
                    axs[2,0].plot(ss[:,49,mi,0],m,markersize=1.)
    for mi,m in enumerate(['ro','yo','ro','yo','go']):
                    axs[0,1].plot(ss[:,0,mi,1],m,markersize=1.)
                    axs[1,1].plot(ss[:,inter,mi,1],m,markersize=1.)
                    axs[2,1].plot(ss[:,49,mi,1],m,markersize=1.)
    for mi,m in enumerate(['bo','yo','bo','yo','go']):
                    axs[0,2].plot(ss[:,0,mi,2],m,markersize=1.)
                    axs[1,2].plot(ss[:,inter,mi,2],m,markersize=1.)
                    axs[2,2].plot(ss[:,49,mi,2],m,markersize=1.)
    axs[0,0].set_ylim(ymin=-0.2, ymax=1.2) 
    axs[1,0].set_ylim(ymin=-0.2, ymax=1.2)    
    axs[2,0].set_ylim(ymin=-0.2, ymax=1.2)    
    axs[0,1].set_ylim(ymin=-0.2, ymax=1.2)    
    axs[1,1].set_ylim(ymin=-0.2, ymax=1.2)    
    axs[2,1].set_ylim(ymin=-0.2, ymax=1.2)    
   

    plt.savefig(ID + "_"+'BIFU.png', bbox_inches='tight',dpi=300, width=40, height=60)
    plt.show()

def bifu_plotG(par,ID,inter):
    IPTG=np.logspace(-2,1,50)
    AHLi=np.logspace(-2,0,100)
    ss = findss(AHLi,IPTG,par)
    fig, axs = plt.subplots(3, 4)
    for mi,m in enumerate(['go','yo','go','yo','ro']):
                    axs[0,0].plot(ss[0,:,mi,0],m,markersize=1.)
                    axs[1,0].plot(ss[inter,:,mi,0],m,markersize=1.)
                    axs[2,0].plot(ss[99,:,mi,0],m,markersize=1.)
    for mi,m in enumerate(['ro','yo','ro','yo','go']):
                    axs[0,1].plot(ss[0,:,mi,1],m,markersize=1.)
                    axs[1,1].plot(ss[inter,:,mi,1],m,markersize=1.)
                    axs[2,1].plot(ss[99,:,mi,1],m,markersize=1.)
    for mi,m in enumerate(['ro','yo','ro','yo','go']):  
                    axs[0,2].plot(ss[0,:,mi,3],m,markersize=1.)
                    axs[1,2].plot(ss[inter,:,mi,3],m,markersize=1.)
                    axs[2,2].plot(ss[99,:,mi,3],m,markersize=1.)
    for mi,m in enumerate(['bo','yo','bo','yo','go']):
                    axs[0,3].plot(ss[0,:,mi,2],m,markersize=1.)
                    axs[1,3].plot(ss[inter,:,mi,2],m,markersize=1.)
                    axs[2,3].plot(ss[99,:,mi,2],m,markersize=1.)

    plt.savefig(ID + "_"+'BIFU.png', bbox_inches='tight',dpi=300, width=40, height=60)
    plt.show()

def bifu_heatmap(p,A,I,ID):
    ss=findss(A,I,p)
    
    hyst_matrix = np.count_nonzero(~np.isnan(ss[:,:,:,0]),axis=2)
    col_matrix = np.nanmax(ss[:,:,:,:],axis=2)
   # hyst_matrix[hyst_matrix==1] = np.NaN
    hyst_matrix = hyst_matrix.astype("float")
    col_matrix = col_matrix.astype("float")
   # col_matrix[col_matrix < .05] = np.NaN
    #hyst_matrix[hyst_matrix==1] = np.NaN
   # hyst_matrix[hyst_matrix==1] = np.NaN
    plt.subplot(1,3,1)
    sns.heatmap(col_matrix[:,:,1], cmap='Reds')
    plt.subplot(1,3,2)
    sns.heatmap(col_matrix[:,:,0], cmap='Greens')
    plt.subplot(1,3,3)
    sns.heatmap(hyst_matrix, cmap='Blues')

    plt.savefig(ID + "_"+'BIFU_heatmap.png', bbox_inches='tight',dpi=300, width=60, height=40)
    plt.show()

def check_par(par1,par2):
#par1='K_GREEN'
#par2='K_RED'
    IPTG=[np.logspace(-2,1,50)[0]]
    AHLi=np.logspace(-2,0,50)

    fig, axs = plt.subplots(10, 10)
    for i,s1 in enumerate(np.arange(2.,step=0.2)+0.2):
        for j,s2 in enumerate(np.arange(2.,step=0.2)+0.2):
            par[par1]=s1
            par[par2]=s2
            ss = findss(AHLi,IPTG,par)
            for mi,m in enumerate(['go','bo','go','ro','yo']):
                axs[i,j].plot(ss[:,0,mi,0],m,markersize=1.)
            axs[i,j].set_ylim(ymin=-0.2,ymax=1.2)
            
    #plt.savefig("plot_TSXL_"+par1+'_' +par2 +'.png', bbox_inches='tight',dpi=300, width=40, height=40)
    plt.show()

def diffusion_plot(u,G,R,par,ID,inter) :
    RED_Giaci = R
    GREEN,RED,AHL,RED_Giac= Integration(G,R,u,RED_Giaci,IPTG,par)
    fig, axs = plt.subplots(3, 2)


    axs[0,0].plot(AHL[0,:,0],'-b')
    axs[0,0].plot(GREEN[0,:,0],'-g')
    axs[0,0].plot(RED[0,:,0],'-r')

    axs[0,1].plot(AHL[-10,:,0],'-b')
    axs[0,1].plot(GREEN[-10,:,0],'-g')
    axs[0,1].plot(RED[-10,:,0],'-r')

    axs[1,0].plot(AHL[0,:,inter],'-b')
    axs[1,0].plot(GREEN[0,:,inter],'-g')
    axs[1,0].plot(RED[0,:,inter],'-r')

    axs[1,1].plot(AHL[-10,:,inter],'-b')
    axs[1,1].plot(GREEN[-10,:,inter],'-g')
    axs[1,1].plot(RED[-10,:,inter],'-r')

    axs[2,0].plot(AHL[0,:,49],'-b')
    axs[2,0].plot(GREEN[0,:,49],'-g')
    axs[2,0].plot(RED[0,:,49],'-r')

    axs[2,1].plot(AHL[-10,:,49],'-b')
    axs[2,1].plot(GREEN[-10,:,49],'-g')
    axs[2,1].plot(RED[-10,:,49],'-r')

    plt.savefig(ID + "_"+'Diffusion.png', bbox_inches='tight',dpi=300, width=60, height=40)



   # for t in np.arange(0,len(AHL[:,0,0])):
    #    axs[0,1].plot(AHL[t,:,24])
       # axs[0,1].set_yscale("log")
    plt.show()


def diffusion_heatmap(p,A,G,R,I,ID,inter,tt):
    RED_Giaci = R
    GREEN,RED,AHL,RED_Giac= Integration(G,R,u,RED_Giaci,IPTG,par,totaltime=tt)
    plt.subplot(3,3,1)
    sns.heatmap(RED[:-10,:,0], cmap='Reds')
    plt.subplot(3,3,2)
    sns.heatmap(GREEN[:-10,:,0], cmap='Greens')
    plt.subplot(3,3,3)
    sns.heatmap(AHL[:-10,:,0], cmap='Blues')
    plt.subplot(3,3,4)
    sns.heatmap(RED[:-10,:,inter], cmap='Reds')
    plt.subplot(3,3,5)
    sns.heatmap(GREEN[:-10,:,inter], cmap='Greens')
    plt.subplot(3,3,6)
    sns.heatmap(AHL[:-10,:,inter], cmap='Blues')
    plt.subplot(3,3,7)
    sns.heatmap(RED[:-10,:,49], cmap='Reds')
    plt.subplot(3,3,8)
    sns.heatmap(GREEN[:-10,:,49], cmap='Greens')
    plt.subplot(3,3,9)
    sns.heatmap(AHL[:-10,:,49], cmap='Blues')

    plt.savefig(ID + "_"+'Diffusion_heatmap.png', bbox_inches='tight',dpi=300, width=60, height=40)
    plt.show()

def diffusion_heatmapG(p,A,G,R,I,ID,inter):
    RED_Giaci = R
    GREEN,RED,AHL,RED_Giac= Integration(G,R,u,RED_Giaci,IPTG,par,totaltime=50)
    plt.subplot(3,3,1)
    sns.heatmap(RED_Giac[:-10,:,0], cmap='Reds')
    plt.subplot(3,3,2)
    sns.heatmap(GREEN[:-10,:,0], cmap='Greens')
    plt.subplot(3,3,3)
    sns.heatmap(AHL[:-10,:,0], cmap='Blues')
    plt.subplot(3,3,4)
    sns.heatmap(RED_Giac[:-10,:,inter], cmap='Reds')
    plt.subplot(3,3,5)
    sns.heatmap(GREEN[:-10,:,inter], cmap='Greens')
    plt.subplot(3,3,6)
    sns.heatmap(AHL[:-10,:,inter], cmap='Blues')
    plt.subplot(3,3,7)
    sns.heatmap(RED_Giac[:-10,:,49], cmap='Reds')
    plt.subplot(3,3,8)
    sns.heatmap(GREEN[:-10,:,49], cmap='Greens')
    plt.subplot(3,3,9)
    sns.heatmap(AHL[:-10,:,49], cmap='Blues')

    plt.savefig(ID + "_"+'Diffusion_heatmap.png', bbox_inches='tight',dpi=300, width=60, height=40)
    plt.show()
#####################################333

IPTG=np.logspace(-2,1,50)

AHLi=np.logspace(-2,0,100)
GREENi=np.ones((100,1))
REDi=np.ones((100,1))
RED_Giaci=np.zeros((100,1))











'''




'''
ID="TSL"
par['K_ahl_green']=20
bifu_plot(par,ID)
bifu_heatmap(par,AHLi,IPTG,ID)


ID="TSLT"
par['K_ahl_green']=1.5
bifu_plot(par,ID)
bifu_heatmap(par,AHLi,IPTG,ID)


ID="TSLT-mushroom"
par['K_ahl_green']=1.5
par['n_ahl_red']=0.4
par['n_ahl_green']=1.4

par['delta_green']=1.2
par['delta_red']=1.4
bifu_plot(par,ID)
bifu_heatmap(par,AHLi,IPTG,ID)


ID="TSXL"
par['K_ahl_green']=20
par['beta_ahl']=1
par['K_ahl_red']=0.4
par['K_ahl']=0.8
inter=10
bifu_plot(par,ID,inter)
bifu_heatmap(par,AHLi,IPTG,ID)
u=np.zeros(25)
G=np.ones((25,1))*0
R=np.ones((25,1))*1
G[10:13]=1
R[10:13]=0
diffusion_heatmap(par,u,G,R,IPTG,ID,inter)


ID="TSXLT"
par['K_ahl_green']=1.5
par['beta_ahl']=1
par['K_ahl']=0.3

par['K_ahl_red']=0.5
par['n_ahl_red']=2
par['n_ahl_green']=2
par['delta_green']=1.
par['delta_red']=1.

#check_par('K_GREEN','K_RED')
inter=16
bifu_plot(par,ID,inter)
bifu_heatmap(par,AHLi,IPTG,ID)
u=np.zeros(30)
G=np.ones((30,1))*1
R=np.ones((30,1))*0
G[10:20]=0
R[10:20]=1
diffusion_heatmap(par,u,G,R,IPTG,ID,inter)



ID="TS-GIAC"
par['K_ahl_green']=20
par['K_ahl_red']=20
par['beta_ahl']=1
par['K_ahl']=1
par['K_GREEN']=1.2
inter=16
bifu_plotG(par,ID,inter)
bifu_heatmap(par,AHLi,IPTG,ID)
u=np.zeros(30)
G=np.ones((30,1))*1
R=np.ones((30,1))*0
G[10:20]=0
R[10:20]=1
diffusion_heatmapG(par,u,G,R,IPTG,ID,inter)



ID="TSXLT_test"
par['K_ahl_green']=1.5
par['beta_ahl']=1
par['K_ahl']=0.3

par['K_ahl_red']=0.5
par['n_ahl_red']=2
par['n_ahl_green']=2
par['delta_green']=1.
par['delta_red']=1.

check_par('K_GREEN','K_ahl_green')

par['K_GREEN']=1
par['K_ahl_green']=0.6


inter=16
bifu_plot(par,ID,inter)
bifu_heatmap(par,AHLi,IPTG,ID)
u=np.zeros(30)
G=np.ones((30,1))*1
R=np.ones((30,1))*0
G[10:20]=0
R[10:20]=1
diffusion_heatmap(par,u,G,R,IPTG,ID,inter)


ID="TSXLT_test2"
par['K_ahl_green']=0.2
par['beta_ahl']=1
par['K_ahl']=1

par['K_ahl_red']=0.2

par['n_ahl_red']=2
par['n_ahl_green']=2
par['delta_green']=1.
par['delta_red']=1.

par['K_GREEN']=0.8
par['K_RED']=0.8

#check_par('n_ahl_red','n_ahl_green')


ss=findss(AHLi,IPTG,par)

inter=16
#bifu_plot(par,ID,inter)
#bifu_heatmap(par,AHLi,IPTG,ID)
u=np.ones(30)*0 #(ss[0,0,4,2]+10e-15)
G=np.ones((30,1))*1#(ss[0,0,4,0]+10e-15)
R=np.ones((30,1))*0#(ss[0,0,4,1]+10e-15)
G[10:20]=0
R[10:20]=1
diffusion_heatmap(par,u,G,R,IPTG,ID,inter)


ID="TSXLT-turing"

par['K_ahl_green']=1.2
par['beta_ahl']=0.5
par['K_ahl']=0
par['K_GREEN']=0.85


par['K_RED']=1.5
par['K_ahl_red']=0
par['n_ahl_red']=2
par['n_ahl_green']=2
par['delta_green']=1.
par['delta_red']=1.
par['delta_ahl']=0.15#0.2
par['D_ahl']=2

inter=10
#bifu_plot(par,ID,inter)
#bifu_heatmap(par,AHLi,IPTG,ID)
u=np.ones(150)*0 #(ss[0,0,4,2]+10e-15)
G=np.ones((150,1))*0#(ss[0,0,4,0]+10e-15)
R=np.ones((150,1))*0#(ss[0,0,4,1]+10e-15)
G[30]=1
diffusion_heatmap(par,u,G,R,IPTG,ID,inter,tt=200)



ID="TSXLT-mushroom2"

par['K_ahl_green']=1.2
par['beta_ahl']=0.5
par['K_ahl']=0
par['K_GREEN']=0.85


par['K_RED']=1.5
par['K_ahl_red']=0
par['n_ahl_red']=2
par['n_ahl_green']=2
par['delta_green']=1.
par['delta_red']=1.
par['delta_ahl']=2.15#0.2
par['D_ahl']=2


inter=10
bifu_plot(par,ID,inter)
bifu_heatmap(par,AHLi,IPTG,ID)
u=np.ones(50)*0 #(ss[0,0,4,2]+10e-15)
G=np.ones((50,1))*0#(ss[0,0,4,0]+10e-15)
R=np.ones((50,1))*0#(ss[0,0,4,1]+10e-15)
G[30:40]=1
diffusion_heatmap(par,u,G,R,IPTG,ID,inter,tt=200)
'''


'''
ID="TSXLT-turing2"

par['K_ahl_green']=1.2
par['beta_ahl']=0.5
par['K_ahl']=0
par['K_GREEN']=0.85


par['K_RED']=1.2
par['K_ahl_red']=0
par['n_ahl_red']=2
par['n_ahl_green']=2
par['delta_green']=1.
par['delta_red']=1.
par['delta_ahl']=0.2#0.2
par['D_ahl']=2


inter=10
bifu_plot(par,ID,inter)
#bifu_heatmap(par,AHLi,IPTG,ID)
u=np.ones(50)*0 #(ss[0,0,4,2]+10e-15)
G=np.ones((50,1))*0#(ss[0,0,4,0]+10e-15)
R=np.ones((50,1))*0#(ss[0,0,4,1]+10e-15)
G[30:40]=1
diffusion_heatmap(par,u,G,R,IPTG,ID,inter,tt=200)   

'''