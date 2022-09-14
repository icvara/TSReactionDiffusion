
## equation for TSLT

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema
import numpy as np
from scipy import optimize
from scipy.optimize import brentq


dtt=0.1
tt=120 #totaltime
#datafile="data_percent.txt"


parlist = [

    {'name':'alpha_red', 'lower_limit':0.,'upper_limit':3.},
    {'name':'basal_red', 'lower_limit':2.,'upper_limit':3.},
    {'name':'beta_red', 'lower_limit':1.,'upper_limit':4.},
    {'name':'K_RED', 'lower_limit':-1.0,'upper_limit':4.0},
    {'name':'n_RED', 'lower_limit':0.0,'upper_limit':2.0},
    {'name':'K_ahl_red', 'lower_limit':-2.0,'upper_limit':4.0},
    {'name':'n_ahl_red', 'lower_limit':0.0,'upper_limit':2.0},
    
    
    {'name':'K_GREEN2', 'lower_limit':-5.0,'upper_limit':3.0},

    {'name':'alpha_green', 'lower_limit':1.,'upper_limit':3.},
    {'name':'basal_green', 'lower_limit':2.,'upper_limit':3.},
    {'name':'beta_green', 'lower_limit':1.,'upper_limit':4.},
    {'name':'K_GREEN', 'lower_limit':0.0,'upper_limit':5.0},
    {'name':'n_GREEN', 'lower_limit':0.,'upper_limit':2.0},
    {'name':'K_ahl_green', 'lower_limit':0.0,'upper_limit':5.0},
    {'name':'n_ahl_green', 'lower_limit':.0,'upper_limit':2.0},
    #{'name':'F_green', 'lower_limit':-4.0,'upper_limit':2.0},
    {'name':'K_IPTG', 'lower_limit':0.0,'upper_limit':5.0},
    #{'name':'K_IPTG2', 'lower_limit':3.0,'upper_limit':8.0}#,
    
    
    {'name':'beta_ahl', 'lower_limit':-1.0,'upper_limit':1.0},
    {'name':'K_ahl', 'lower_limit':-4.0,'upper_limit':1.0},
    {'name':'n_ahl', 'lower_limit':0.0,'upper_limit':2.0}
  
    
]






#first part on simualtion
########################################################

def model_TSXLT(GREENi,REDi,REDfi,AHLi,A0,IPTG,par):

    #GREENi=np.maximum(GREENi - 10**par['basal_green'],0)
    #REDi=np.maximum(REDi - 10**par['basal_red'],0)

    AHL_tot = AHLi + A0#[:,np.newaxis,np.newaxis]

    GREEN = 10**par['alpha_green'] + (10**par['beta_green']*np.power(AHL_tot*10**par['K_ahl_green'],par['n_ahl_green']))/(1+np.power(AHL_tot*10**par['K_ahl_green'],par['n_ahl_green']))
    GREEN = GREEN / (1 + np.power(REDi*10**par['K_RED'],par['n_RED']))
    GREEN = GREEN - par['delta_green']*GREENi 

    free_GREENi= GREENi / ( 1+ 10**par['K_IPTG']*IPTG)#[np.newaxis,:,np.newaxis]

    RED = 10**par['alpha_red'] + (10**par['beta_red']*np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))
    RED = RED / (1 + np.power(free_GREENi*10**par['K_GREEN'],par['n_GREEN']))
    RED = RED - par['delta_red']*REDi # + par['alpha_red']

    REDf = 10**par['alpha_red'] + (10**par['beta_red']*np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))
    REDf = REDf / (1 + np.power(free_GREENi*10**par['K_GREEN2'],par['n_GREEN']))
    REDf = REDf - par['delta_red']*REDfi # + par['alpha_red']

    AHL = (10**par['beta_ahl']*np.power(GREENi*10**par['K_ahl'],par['n_ahl']))/(1+np.power(GREENi*10**par['K_ahl'],par['n_ahl']))
    AHL = AHL - par['delta_ahl']*(AHLi) 

    #GREEN=GREEN + 10**par['basal_green']
    #REDf=REDf + 10**par['basal_red']

    return GREEN,RED,REDf,AHL



################################3
### dynamic analysis
########################################


   
def solvedfunction(Gi,A,I,par,model):
    #rewrite the system equation to have only one unknow and to be call with scipy.optimze.brentq
    #the output give a function where when the line reach 0 are a steady states
    par['delta_green']=1
    par['delta_red']=1 #1.2
    par['delta_ahl']=1 #1.2


    if model == 'TSLT':

        Gf = Gi / ( 1+ 10**par['K_IPTG']*I)

        R = 10**par['alpha_red'] + ( 10**par['beta_red']*np.power(A*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(A*10**par['K_ahl_red'],par['n_ahl_red']))
        R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN']))  
        R = ( R ) / par['delta_red']  


        G = 10**par['alpha_green'] + ( 10**par['beta_green']*np.power(A*10**par['K_ahl_green'],par['n_ahl_green']))/(1+np.power(A*10**par['K_ahl_green'],par['n_ahl_green'])) 
        G = G / (1 + np.power(R*10**par['K_RED'],par['n_RED'])) 
        G = (G ) / par['delta_green']

        func = G - Gi


    elif model == 'TSXLT':
      
      
        Gf = Gi / ( 1+ 10**par['K_IPTG']*I)
        
        AHL =  ( 10**par['beta_ahl']*np.power(Gi*10**par['K_ahl'],par['n_ahl']))/(1+np.power(Gi*10**par['K_ahl'],par['n_ahl'])) 
        AHL = (AHL) / par['delta_ahl']
        AHL_tot = A + AHL

        R = 10**par['alpha_red'] + ( 10**par['beta_red']*np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))
        R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN']))  
        R = ( R ) / par['delta_red']  


        G = 10**par['alpha_green'] + ( 10**par['beta_green']*np.power(AHL_tot*10**par['K_ahl_green'],par['n_ahl_green']))/(1+np.power(AHL_tot*10**par['K_ahl_green'],par['n_ahl_green']))
        G = G / (1 + np.power(R*10**par['K_RED'],par['n_RED'])) 
        G = (G ) / par['delta_green']

        func = G - Gi

    else:
        print("wrong model type, try TSLT or TSXLT. you type: " +model)
        func=0

    return func 


def findss(A,I,par,model):
    #list of fixed par
    #function to find steady state
    #1. find where line reached 0
    ss=[]
    nNode=3 # number of nodes : X,Y,
    if model == 'TSXLT':
        nNode = 4
    nStstate= 7
    nAHL= len(A)
    nIPTG=len(I)
    ss=np.ones((nAHL,nIPTG,nStstate,nNode))*np.nan  
    for ai,a in enumerate(A):
        for iptgi,iptg in enumerate(I):
            Gi=np.logspace(-50,10,10000,base=10)
            f=solvedfunction(Gi,a,iptg,par,model)         
            x=f[1:-1]*f[0:-2] #when the output give <0, where is a change in sign, meaning 0 is crossed
            index=np.where(x<0)
            for it,i in enumerate(index[0]):
                if model == 'TSLT':
                    G = brentq(solvedfunction, Gi[i], Gi[i+1],args=(a,iptg,par,model)) #find the value of G at 0
                    Gf = G / ( 1+ 10**par['K_IPTG']*iptg)
                    R = 10**par['alpha_red'] + ( 10**par['beta_red']*np.power(a*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(a*10**par['K_ahl_red'],par['n_ahl_red'])) 
                    #K_GREEN change to F_GREEN, try to change dynamic from TETR and mcherry repression by LacI
                    #R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN']))
                    R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN'])) 
                    R = ( R ) / par['delta_red'] 

                    Rf = 10**par['alpha_red'] + ( 10**par['beta_red']*np.power(a*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(a*10**par['K_ahl_red'],par['n_ahl_red'])) 
                    Rf = Rf / (1 + np.power(Gf*10**par['K_GREEN2'],par['n_GREEN'])) 
                    Rf = ( Rf ) / par['delta_red'] 


                    Rf = Rf + 10**par['basal_red'] #autofluo
                    G = G + 10**par['basal_green'] #autofluo
                    ss[ai,iptgi,it]=np.array([G,Rf,R])


                elif model == 'TSXLT':
                    
                    G = brentq(solvedfunction, Gi[i], Gi[i+1],args=(a,iptg,par,model)) #find the value of G at 0
                    Gf = G / ( 1+ 10**par['K_IPTG']*iptg)
                   

                    AHL =  ( 10**par['beta_ahl']*np.power(G*10**par['K_ahl'],par['n_ahl']))/(1+np.power(G*10**par['K_ahl'],par['n_ahl'])) 
                    AHL = (AHL) / par['delta_ahl']
                    AHL_tot = a + AHL

                    R = 10**par['alpha_red'] + ( 10**par['beta_red']*np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red'])) 
                    R = R / (1 + np.power(Gf*10**par['K_GREEN'],par['n_GREEN'])) 
                    R = ( R ) / par['delta_red'] 


                    Rf = 10**par['alpha_red'] + ( 10**par['beta_red']*np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red']))/(1+np.power(AHL_tot*10**par['K_ahl_red'],par['n_ahl_red'])) 
                    Rf = Rf / (1 + np.power(Gf*10**par['K_GREEN2'],par['n_GREEN'])) 
                    Rf = ( Rf ) / par['delta_red'] 


                    Rf = Rf + 10**par['basal_red'] #dissociate TetR form mcherry
                    G = G + 10**par['basal_green'] #dissociate LacI form GFP
                    
                    ss[ai,iptgi,it]=np.array([G,Rf,R,AHL])

             #   ss[ai,iptgi,it]=np.array([G+par['basal_green'],R+par['basal_red']])
    return ss




def jacobianMatrix(ss,par,A0,I0,model="TSXLT"):
    #doesnt work... not goood ...
    print("use approximateJacob")

    '''
    JM=np.ones((ss.shape[0],ss.shape[1],ss.shape[2],3,3))*np.nan 
    if model=="TSXLT":
        G=np.copy(ss[:,:,:,0]) - 10**par['basal_green']
        #Rf=np.copy(ss[:,:,:,1])
        R=np.copy(ss[:,:,:,2])
        A=np.copy(ss[:,:,:,3])
        I= I0[np.newaxis,:,np.newaxis]

        A=A+A0[:,np.newaxis,np.newaxis]

        Gf = G / ( 1+ 10**par['K_IPTG']*I0[np.newaxis,:,np.newaxis])

        #IPTG somewheere...

        dGdg = - par['delta_green'] *G/G

        dGdr = -(((10**par['beta_green']*np.power((10**par['K_ahl_green']*A),par['n_ahl_green']))/(1+np.power((10**par['K_ahl_green']*A),par['n_ahl_green']))+10**par['alpha_green'])*par['n_RED']*np.power((10**par['K_RED']*R),par['n_RED']))
        dGdr =  dGdr/(R*np.power((np.power((10**par['K_RED']*R),par['n_RED'])+1),2))

        dGda = 10**par['beta_green']*par['n_ahl_green']*np.power((A*10**par['K_ahl_green']),par['n_ahl_green'])
        dGda= dGda/ ((np.power(R*10**par['K_RED'],par['n_RED'])+1)*A*np.power((np.power(A*10**par['K_ahl_green'],par['n_ahl_green'])+1),2))

        dRfdr = - par['delta_red'] *Rf/Rf
        
        dRfdg = -(((par['beta_red']*np.power((10**par['K_ahl_red']*A),par['n_ahl_red']))/(1+np.power((10**par['K_ahl_red']*A),par['n_ahl_red']))+par['alpha_red'])*par['n_GREEN']*np.power((10**par['K_GREEN2']*Gf),par['n_GREEN']))
        dRfdg =  dRfdg/(Gf*np.power((np.power((10**par['K_GREEN2']*Gf),par['n_GREEN'])+1),2))

        dRfda = par['beta_red']*par['n_ahl_red']*np.power((A*10**par['K_ahl_red']),par['n_ahl_red'])
        dRfda= dRfda/ ((np.power(G*10**par['K_GREEN2'],par['n_GREEN'])+1)*A*np.power((np.power(A*10**par['K_ahl_red'],par['n_ahl_red'])+1),2))
 
        dRdr = - par['delta_red'] *R/R


       # dRdg = - (10**par['alpha_red']+10**par['beta_red'])*np.power((10**par['K_ahl_red']*A),par['n_ahl_red'])*par['n_GREEN']*np.power((10**par['K_GREEN']*G/(10**par['K_IPTG']*I+1)),par['n_GREEN'])
       # dRdg = dRdg / (np.power((10**par['K_ahl_red']*A),par['n_ahl_red'])+1)*G*np.power((np.power((10**par['K_GREEN']*G/(10**par['K_IPTG']*I+1)),par['n_GREEN'])+1),2)
        dRdg = -(((10**par['beta_red']*np.power((10**par['K_ahl_red']*A),par['n_ahl_red']))/(1+np.power((10**par['K_ahl_red']*A),par['n_ahl_red']))+10**par['alpha_red'])*par['n_GREEN']*np.power((10**par['K_GREEN']*Gf),par['n_GREEN']))
        dRdg =  dRdg/(Gf*np.power((np.power((10**par['K_GREEN']*Gf),par['n_GREEN'])+1),2))


       # dRda= - (10**par['alpha_red']-10**par['beta_red'])*par['n_ahl_red']*np.power((10**par['K_ahl_red']*A),par['n_ahl_red'])
       # dRda= dRda / (np.power((10**par['K_GREEN']*G/(10**par['K_IPTG']*I+1)),par['n_GREEN']+1)*A*np.power((np.power((10**par['K_ahl_red']*A),par['n_ahl_red'])+1),2))

        dRda = 10**par['beta_red']*par['n_ahl_red']*np.power((A*10**par['K_ahl_red']),par['n_ahl_red'])
        dRda= dRda/ ((np.power(Gf*10**par['K_GREEN'],par['n_GREEN'])+1)*A*np.power((np.power(A*10**par['K_ahl_red'],par['n_ahl_red'])+1),2))

        dAdg = 10**par['beta_ahl']*par['n_ahl']*np.power((G*10**par['K_ahl']),par['n_ahl'])
        dAdg = dAdg/ (G*np.power((np.power((10**par['K_ahl']*G),par['n_ahl'])+1),2))

        dAdr = 0 * A/A

        dAda =  - par['delta_ahl'] *A/A
           
        JM[:,:,:,0,0]=dGdg
        JM[:,:,:,0,1]=dGdr
        JM[:,:,:,0,2]=dGda

        JM[:,:,:,1,0]=dRfdg
        JM[:,:,:,1,1]=dRfdr
        JM[:,:,:,1,2]=dRfda

        JM[:,:,:,1,0]=dRdg
        JM[:,:,:,1,1]=dRdr
        JM[:,:,:,1,2]=dRda

        JM[:,:,:,2,0]=dAdg
        JM[:,:,:,2,1]=dAdr
        JM[:,:,:,2,2]=dAda
    else:
        JM=0   
    return JM
    '''
    return 0


def approximateJacob(G,R,Rf,A,A0,I,par): 
    #allow to check if jacobian matrix derivate are correctly written
    delta=10e-5
    g,rf,r,a = model_TSXLT(G,R,Rf,A,A0,I,par)


    dgdg = (model_TSXLT(G+delta,R,Rf,A,A0,I,par)[0]-g)/delta
    dgdr = (model_TSXLT(G,R+delta,Rf,A,A0,I,par)[0]-g)/delta
    dgda = (model_TSXLT(G,R,Rf,A+delta,A0,I,par)[0]-g)/delta

    drdg = (model_TSXLT(G+delta,R,Rf,A,A0,I,par)[1]-r)/delta
    drdr = (model_TSXLT(G,R+delta,Rf,A,A0,I,par)[1]-r)/delta
    drda = (model_TSXLT(G,R,Rf,A+delta,A0,I,par)[1]-r)/delta

    dadg = (model_TSXLT(G+delta,R,Rf,A,A0,I,par)[3]-a)/delta
    dadr = (model_TSXLT(G,R+delta,Rf,A,A0,I,par)[3]-a)/delta
    dada = (model_TSXLT(G,R,Rf,A+delta,A0,I,par)[3]-a)/delta

    #JM=np.ones((G.shape[0],3,3))*np.nan 
    JM=np.ones((G.shape[0],G.shape[1],G.shape[2],3,3))*np.nan 


    JM[:,:,:,0,0]=dgdg
    JM[:,:,:,0,1]=dgdr
    JM[:,:,:,0,2]=dgda

    JM[:,:,:,1,0]=drdg
    JM[:,:,:,1,1]=drdr
    JM[:,:,:,1,2]=drda

    JM[:,:,:,2,0]=dadg
    JM[:,:,:,2,1]=dadr
    JM[:,:,:,2,2]=dada

    return JM



#findss(A,I,par,model):

def getEigen(ss,A0,I0,par,model="TSXLT"):
    if model !="TSXLT":
        return 0
    #JA=jacobianMatrix(ss,par,A0,I0,model)
    A0=A0[:,np.newaxis,np.newaxis]
    I0=I0[np.newaxis,:, np.newaxis]
    J=approximateJacob(ss[:,:,:,0]-10**par['basal_green'],ss[:,:,:,2],ss[:,:,:,1]-10**par['basal_red'],ss[:,:,:,3],A0,I0,par)

   # print(J-JA)
    J2=np.nan_to_num(J)
    eigvals, eigvecs =np.linalg.eig(J2)
    sse=eigvals
    return sse #, np.trace(A), np.linalg.det(A)



def TuringInstability(A0,I0,par,n=100,model="TSXLT"):
    if model !="TSXLT":
        return 0
    q=np.linspace(0,20,n)#100


    ss=findss(A0,I0,par,model)
    #J=jacobianMatrix(ss,par,A0,model)
    A=A0[:,np.newaxis,np.newaxis]
    I=I0[np.newaxis,:,np.newaxis]
    J=approximateJacob(ss[:,:,:,0]-10**par['basal_green'],ss[:,:,:,2],ss[:,:,:,1]-10**par['basal_red'],ss[:,:,:,3],A,I,par)

    JE=np.ones((len(A0),len(I0),ss.shape[2],n,3,3))*np.nan
    JE[:,:,:,:,:,:]=np.copy(J[:,:,:,np.newaxis,:,:])
    #JE=np.ones((len(A0),len([I]),5,3,3,20))*np.copy(J)

    eigens=np.ones((len(A0),len([I]),5,len(q),3))*np.nan


    if model =="TSXLT":
        JE[:,:,:,:,2,2]=JE[:,:,:,:,2,2]- q**2*par['D_ahl']
        JE2=np.nan_to_num(JE)
        eigvals, eigvecs =np.linalg.eig(JE2)
        sse=eigvals
        sse=np.apply_along_axis(sortComplex,4,sse)


    return sse

def sortComplex(a):
    b=np.sort_complex(a)


    return b



###########################################



def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def Get_data(dataname, st='median'):
    path=dataname
    df = pd.read_csv(path,sep='\t' ,header=[0])
    df[df == ' NA'] = np.nan

    t=find_between(dataname, '_','.')
    #t=find_between(dataname, '_','_')
   

    if t=="gated":
        names=df.columns
        df_G=df[df['fluo'] == "GREEN"]
        df_gg=df_G[df_G["minmax"] == "max"]
        df_gr=df_G[df_G["minmax"] == "min"]
        gg = df_gg.pivot(index='AHL', columns='IPTG', values='mean').astype(float)
        gr = df_gr.pivot(index='AHL', columns='IPTG', values='mean').astype(float)

        df_R=df[df['fluo'] == "RED"]
        df_rg=df_R[df_R["minmax"] == "min"]
        df_rr=df_R[df_R["minmax"] == "max"]
        rg = df_rg.pivot(index='AHL', columns='IPTG', values='mean').astype(float)
        rr = df_rr.pivot(index='AHL', columns='IPTG', values='mean').astype(float)

    elif t=="percent":
        names=df.columns
        df_G=df[df['gate'] == 1]
        df_gg=df_G[df_G["sample"] == "G"]
        df_gr=df_G[df_G["sample"] == "R"]
        gg = df_gg.pivot(index='AHL', columns='IPTG', values='mean').astype(float)
        gr = df_gr.pivot(index='AHL', columns='IPTG', values='mean').astype(float)

        df_R=df[df['gate'] == 2]
        df_rg=df_R[df_R["sample"] == "G"]
        df_rr=df_R[df_R["sample"] == "R"]
        rg = df_rg.pivot(index='AHL', columns='IPTG', values='mean').astype(float)
        rr = df_rr.pivot(index='AHL', columns='IPTG', values='mean').astype(float)

    #elif t=="":
    else:
        names=df.columns
        df_G=df[df['fluo'] == "GREEN"]
        df_gg=df_G[df_G["sample"] == "G"]
        df_gr=df_G[df_G["sample"] == "R"]
        gg = df_gg.pivot(index='AHL', columns='IPTG', values=st).astype(float)
        gr = df_gr.pivot(index='AHL', columns='IPTG', values=st).astype(float)

        df_R=df[df['fluo'] == "RED"]
        df_rg=df_R[df_R["sample"] == "G"]
        df_rr=df_R[df_R["sample"] == "R"]
        rg = df_rg.pivot(index='AHL', columns='IPTG', values=st).astype(float)
        rr = df_rr.pivot(index='AHL', columns='IPTG', values=st).astype(float)

    return  gg,gr,rr,rg # gg = green fluo, green sample


###########################################3


def distance(pars,path,modeltype):
    
   # GG,GR,GA,RG,RR,RA = model(pars,totaltime, dt)
    gg,gr,rr,rg=Get_data(path)
    AHL=gg.index.values
    IPTG=gg.columns.values
    
    pars['delta_green']=1
    pars['delta_red']=1
    ss= findss(AHL,IPTG,pars,modeltype)

    M=np.nanmax(ss[:,:,:,:],axis=2)
    m=np.nanmin(ss[:,:,:,:],axis=2)
    
    
    d_green = np.sqrt(np.power(gg.to_numpy() - M[:,:,0],2))
    d_red = np.sqrt(np.power(rr.to_numpy() - M[:,:,1],2))

    d_green2 = np.sqrt(np.power(gr.to_numpy() - m[:,:,0],2))
    d_red2 = np.sqrt(np.power(rg.to_numpy() - m[:,:,1],2))
    
    d=np.nansum((d_green,d_red,d_red2,d_green2))
    #d=np.nansum((d_green,d_green2))
    #d=np.nansum((d_red,d_red2))

    return d


#############################################3



