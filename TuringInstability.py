'''
===================================

Script to find Turing instability

================================
'''

import numpy as np
import time
import multiprocessing
from functools import partial
import random
from model import *



'''
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy.linalg import eig
#from scipy.stats import norm, uniform, multivariate_normal
from scipy import optimize
from scipy.optimize import brentq
import pandas as pd
'''

initdist=1000000000
finaldist=100.

version="MAP_TuringInstability"

if os.path.isdir(version) is False: ## if 'smc' folder does not exist:
        os.mkdir(version) ## create it, the output will go there

path='/users/ibarbier/RD/'
#path='C:/Users/Administrator/Desktop/Modeling/RD/'

parlist = [ 
    #GREEN
    {'name' : 'K_ahl_green', 'lower_limit':-1.0,'upper_limit':2.0},
    {'name' : 'n_ahl_green', 'lower_limit':0.5,'upper_limit':2.0},
    {'name' : 'alpha_green', 'lower_limit':0.001,'upper_limit':1.0},
    {'name' : 'beta_green', 'lower_limit':0.0,'upper_limit':4.0},
 #   {'name' : 'delta_green', 'lower_limit':1.0,'upper_limit':1.0},
    {'name' : 'K_GREEN', 'lower_limit':0.0,'upper_limit':4.0},
   # {'name' : 'n_GREEN', 'lower_limit':0.5,'upper_limit':2.0},

    #RED
    {'name' : 'K_ahl_red', 'lower_limit':-1.0,'upper_limit':3.0},
    {'name' : 'n_ahl_red', 'lower_limit':0.5,'upper_limit':2.0},
    {'name' : 'alpha_red', 'lower_limit':0.0,'upper_limit':1.0},
    {'name' : 'beta_red', 'lower_limit':0.0,'upper_limit':4.0},
 #   {'name' : 'delta_red', 'lower_limit':1.0,'upper_limit':1.0},
    {'name' : 'K_RED', 'lower_limit':0.0,'upper_limit':4.0},
  #  {'name' : 'n_RED', 'lower_limit':0.5,'upper_limit':2.0},

    #AHL
    {'name' : 'K_ahl', 'lower_limit':0.0,'upper_limit':1.0},
    {'name' : 'n_ahl', 'lower_limit':0.5,'upper_limit':2.0},
    {'name' : 'beta_ahl', 'lower_limit':0.0,'upper_limit':1.0},
  #  {'name' : 'delta_ahl', 'lower_limit':1.0,'upper_limit':1.0},
  #  {'name' : 'D_ahl', 'lower_limit':0.01,'upper_limit':1.0},
]

def addfixedpar(par):
    #list of fixed par
    par['delta_red']=1
    par['delta_green']=1
    par['delta_ahl']=1
    par['D_ahl']=1
    par['K_IPTG']=0
    par['n_GREEN']=2
    par['n_RED']=2

    return par



def choosepar(parlist):
    #choose random par in the defined range
    samplepar=[]
    for ipar,par in enumerate(parlist):
        samplepar.append(random.uniform(par['lower_limit'], par['upper_limit']))
   # p=pars_to_dict(samplepar,parlist)
    return np.array(samplepar)


def pars_to_dict(pars,parlist):
### This function is not necessary, but it makes the code a bit easier to read,
### it transforms an array of pars e.g. p[0],p[1],p[2] into a
### named dictionary e.g. p['k0'],p['B'],p['n'],p['x0']
### so it is easier to follow the parameters in the code
    dict_pars = {}
    for ipar,par in enumerate(parlist):
        dict_pars[par['name']] = pars[ipar] 
    return dict_pars






def isOscillation(e):
   b=0
   if np.any(e.real>0): #unstable
               #check complex conjugate
      if len(e.real) != len(set(e.real)):
         if np.any(e.imag==0) and len(e.imag) != len(set(np.abs(e.imag))):
            b=1
   return b

def isStability(e):
   b=0
   if np.all(e.real<0): #stable
            b=1
   return b



def getTuringinstability(par,n=100,model="TSXLT"):
   A=np.zeros(1)

   par=addfixedpar(par)
   turing_type= -1 #no stability ? 
   newpar=[]
   seeD = TuringInstability(A,par,n,model)
   '''
   three case
   1: stable becomes unstable
   2: stable becomes unstable and then stable again
   3: oscillation becomes stable

   '''
   oscil=np.apply_along_axis(isOscillation,3,seeD)
   stab=np.apply_along_axis(isStability,3,seeD)

   if np.any(stab==1):
      indexa=np.where(stab==1)
      a=seeD[0,0,indexa[2][0],:,:] 
      x=a[1:-1,:]*a[0:-2,:] #when the output give <0, where is a change in sign, meaning 0 is crossed
      index=np.where(x<0)
      if len(set(index[0]))==0: #no change in signe
         turing_type=0
      if len(set(index[0]))==1: # change in signe once
         turing_type=1
      if len(set(index[0]))>1: # change in signe at least twice. classic Turing pattern
         turing_type=2

   if np.any(oscil==1):
      indexb=np.where(oscil==1)
      b=seeD[0,0,indexb[2][0],:,:]
      stab_b=np.apply_along_axis(isStability,1,b)
      if np.any(stab_b==1):
         turing_type=3
                  
   return turing_type


def choosepar(parlist):
    #choose random par in the defined range
    samplepar=[]
    for ipar,par in enumerate(parlist):
        samplepar.append(random.uniform(par['lower_limit'], par['upper_limit']))
   # p=pars_to_dict(samplepar,parlist)
    return np.array(samplepar)

def GeneratePars(parlist, ncpus,Npars=1000):
    #@EO: to compare the 2 versions, set the seed
    #np.random.seed(0)

    ## Call to function GeneratePar in parallel until Npars points are accepted
    trials = 0
    start_time = time.time()
    results = []

    pool = multiprocessing.Pool(ncpus)
    results = pool.map(func=partial(calculatePar, parlist),iterable=range(Npars), chunksize=10)
    pool.close()
    pool.join()    
    end_time = time.time()
    print(f'>>>> Loop processing time: {end_time-start_time:.3f} sec on {ncpus} CPU cores.')    

    newparlist = [result[0] for result in results]
    turingtype = [result[1] for result in results]

    return(newparlist,turingtype)


def calculatePar(parlist, iter):
  #selectpar=[]
  newpar=choosepar(parlist)    
  p=pars_to_dict(newpar,parlist)
  tutype = getTuringinstability(p)
  #if tu >0:
    #selectpar.append(newpar)
  return newpar,tutype


def load(name,parlist):
    tutype= np.loadtxt(name+"_turingtype.out")
    p= np.loadtxt(name+"_par.out")

    namelist=[]
    for i,par in enumerate(parlist):
        namelist.append(parlist[i]['name'])
    df = pd.DataFrame(p, columns = namelist)
    df['tutype']=tutype
    df=df.sort_values(by='tutype', ascending=False)
    return df

def run(name,Npars=5000,ncpus=40):
    par,tutype=GeneratePars(parlist, ncpus,Npars=Npars)
    np.savetxt(version +"/"+ name+'_turingtype.out', tutype)
    np.savetxt(version +"/"+name+'_par.out', par)


print(tutype)
run("test",10,ncpus=1)