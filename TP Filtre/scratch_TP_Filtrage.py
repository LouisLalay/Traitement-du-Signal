#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 09:05:45 2019

@author: said
"""

#%% sections de modules à importer 

import numpy as np
import scipy.signal
import soundfile as sf

import matplotlib.pyplot as plt
import os
import tempfile
import platform

if platform.system()=='Darwin': #MAC 
    import sounddevice as sd
    
def play(y,Fe=44100):
    z=np.real(y)/(abs(np.real(y)).max())
    if platform.system()=='Darwin': #MAC (sous linux sounddevice a un comportement erratique)
        sd.play(z,Fe)
        return
    fichier=tempfile.mktemp()+'SON_TP.wav'
    sec=len(y)/Fe
    if sec<=20:
        rep=True
    if sec>20:
        print ('Vous allez créer un fichier son de plus de 20 secondes.')
        rep=None
        while rep is None:
            x=input('Voulez-vous continuer? (o/n)')
            if x=='o':
                rep=True
            if x=='n':
                rep=False
            if rep is None:
                print ('Répondre par o ou n, merci. ')
    if rep:
        sf.write(fichier,z,Fe)
        os.system('/usr/bin/play '+fichier+' &')


#%%
3+9

np.arange(0, 1, 0.1)
np.zeros((2, 2))
np.zeros((3, 1)) # colonne

np.zeros((1, 3)) #ligne

#%%
A=np.asarray([[1,2],[ 3 ,4]])
B=np.asarray([[-2,1],[ 1.5 ,-0.5]])
print('A=\n',A)
print('B=\n',B)
print('produit matriciel\n', A@B)
print('produit point a point\n', A*B)

#%% RESHAPE
A=np.asarray([[1,2],[3,4]])
A.reshape(-1)
u=np.asarray([[1, 2, 2, 1]]) # matrice ligne
v=np.asarray([[1],[2],[3],[4]]) # matrice colonne
print('u*v', u*v) # résultat inattendu (produit tensoriel)
print('u.reshape(-1)*v.reshape(-1)',u.reshape(-1)*v.reshape(-1)) #on est s^ur de la forme c’est le produit point a point

#%%
def valconv(h,u): # cette version ne fonctionne pas telle quelle!
    a=len(h) # pour un tableau numpy len(v)==v.shape[0]
    b=len(u)
    u=u.reshape(-1) #s’assurer que u et v son de la meme forme
    h=h.reshape(-1) #
    out = np.zeros(a+b-1)
   
    idx=None
    
    out[n]=   (h[idx]*u[n-idx]).sum()
    return out



#%% ATTENTION de bien mettre 
repertoire='/home/user/' #mettre le chemin ou vous avez mis les données
[x,Fe]=sf.read(repertoire+'oiseau.wav')

play(x,Fe)
#%%
fichier='handel.wav'
#fichier='oiseau.wav'
#fichier='marche.wav'
[x,Fe]=sf.read(repertoire+fichier)
play(x,Fe)
#%%
h=np.zeros(2*int(Fe/5))
h[-1]=0.64
h[int(Fe/3)]=0.8
if len(x.shape)>1: # prendre un seul canal s'il y en a deux
    y=x[:,0].copy()
else:
    y=x.copy()
yecho=scipy.convolve(y,h)
play(np.concatenate((y,yecho),Fe))

#%%

def rejette1(entree,f0,Fs):

#valeur recommandée pour rho: au moins 0.9 et jamais plus de 1
    lfilter=scipy.signal.lfilter
    z0=exp(2*i*pi*f0/Fs); # position sur le cercle unite de la fréquence à éliminer
    
    sortie=lfilter([1 ,-z0],1,entree); 
    
    sortie=lfilter([1, -np.conj(z0)],1,sortie); #Le signal est reel, il faut donc éliminer -f0 aussi!
    
    
    sortie=np.real(sortie);# En raison d'erreurs d'arrondi il peut subsister une partie imaginaire
    return sortie
def rejette2(entree,f0,Fs,rho):

#valeur recommandée pour rho: au moins 0.9 et jamais plus de 1
    lfilter=scipy.signal.lfilter
    z0=exp(2*i*pi*f0/Fs); # position sur le cercle unite de la fréquence à éliminer
    
    sortie=lfilter([1 ,-z0],[1 ,-rho*z0],entree); 
    
    sortie=lfilter([1, -np.conj(z0)],[1 ,-rho*np.conj(z0)],sortie); #Le signal est reel, il faut donc éliminer -f0 aussi!

    sortie=np.real(sortie);# En raison d'erreurs d'arrondi il peut subsister une partie imaginaire
    return sortie





#%% parasite
y=x[:4*Fe,0]
n=np.arange(len(y)); #le temps discret que dure le signal y
f0=1261;
yb=y+0.1*np.cos(2*pi*f0/Fe*n); #on ajoute une onde parasite
play(yb,Fe)
#%% nettoyage par simple rejection
energ=(y**2).sum()**0.5
ynet=rejette1(yb,f0,Fe)
energn=(ynet**2).sum()**0.5
play(ynet/energn*energ,Fe) #normalisation
#%% nettoyage en ajoutant un pole dans le filtre
energ=(y**2).sum()**0.5
ynet=rejette2(yb,f0,Fe,0.99) # faire évoluer rho
energn=(ynet**2).sum()**0.5
play(ynet/energn*energ,Fe) #normalisation

#%%
def norm2(x):
    return (x**2).sum()**0.5
#%%
def liste_positions_symetriques(Tx=10,Ty=5,Tz=2.5,xsource=3,ysource=2.5,zsource=1,nrebonds=8):

    accu=[]
    listetrav=[(xsource,ysource,zsource)]
    accu=accu+[(xsource,ysource,zsource)]
    for k in range(nrebonds):

        rebx0= [(-m[0]    ,m[1],m[2]) for m in listetrav] #rebonds sur le mur X=0
        rebxTx=[(2*Tx-m[0],m[1],m[2]) for m in listetrav] #rebonds mur X=Tx
        
        reby0= [(m[0],-m[1]    ,m[2]) for m in listetrav] #rebonds sur le mur Y=0
        rebyTy=[(m[0],2*Ty-m[1],m[2]) for m in listetrav] #rebonds mur Y=Ty

        rebz0= [(m[0],m[1],-m[2]    ) for m in listetrav] #rebonds sur le mur Z=0
        rebzTz=[(m[0],m[1],2*Tz-m[2]) for m in listetrav] #rebonds mur Z=Tz

        listetrav=rebx0+rebxTx+reby0+rebyTy
        accu+=listetrav
    return accu
#%% 
def echosmultiples(Tx=10,Ty=5,Tz=2.5,xsource=3,ysource=2.5,zsource=1,xauditeur=8,yauditeur=1,zauditeur=1,nrebonds=8):
    c=300 # metres/seconde
    listepositions=liste_positions_symetriques(Tx=Tx,Ty=Ty,Tz=Tz,xsource=xsource,ysource=ysource,nrebonds=nrebonds)
    d=[((m[0]-xauditeur)**2+(m[1]-yauditeur)**2+(m[1]-zauditeur)**2)**0.5 for m in listepositions]
    h=np.zeros(int(max(d)/c*Fe+1))

    for k in d:
        posdansh=int(np.floor(k/c*Fe))
        h[posdansh]+=1/k
    return h

#%%
h=echosmultiples(Tx=30,Ty=30,Tz=3,xsource=5,ysource=5,zsource=0.5,xauditeur=20,yauditeur=2,zauditeur=1,nrebonds=8) 
#%%
if len(x.shape)>1:
    y=x[:,0].copy()
else:
    y=x.copy()
ny=norm2(y)
yechomult=scipy.convolve(y,h)


play(np.concatenate((y,yechomult)))
