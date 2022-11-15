import numpy as np

def valconv(h,u,n):
    a=len(h) # pour un tableau numpy len(v)==v.shape[0]
    b=len(u)
    u=u.reshape(-1) #s'assurer que u et v son de la meme forme
    h=h.reshape(-1)
    if n == 0:
        return h[0]*u[0]
    idx = np.arange(max(0,n-b+1), min(a-1, n) + 1)
    return (h[idx]*u[n-idx]).sum()

print(valconv(np.ones(3),np.ones(2),0))
print(valconv(np.ones(3),np.ones(2),2))

def full_conv(h, u):
    a=len(h) # pour un tableau numpy len(v)==v.shape[0]
    b=len(u)
    out = np.zeros((a+b-1))
    u=u.reshape(-1) #s'assurer que u et v son de la meme forme
    h=h.reshape(-1)
    out[0] = h[0]*u[0]
    for n in range(1, a+b-1):
        idx = np.arange(max(0,n-b+1), min(a-1, n) + 1)
        out[n] = (h[idx]*u[n-idx]).sum()
    return out

print(full_conv(np.ones(5),np.ones(7)))
print(np.convolve(np.ones(5),np.ones(7)))

import soundfile as sf
import matplotlib.pyplot as plt
import sounddevice as sd
from time import time

[x,Fe]=sf.read('piano.wav')
def echo(x, Fe, t1):
    n1 = round(t1*Fe)
    h = np.zeros((n1+1))
    h[0]=1
    h[n1]=0.8
    return full_conv(h, x)

def echo_scipy(x, Fe, t1):
    n1 = round(t1*Fe)
    h = np.zeros((n1+1))
    h[0]=1
    h[n1]=0.8
    return np.convolve(h, x)
t1 = 0.01
start = time()
y = echo(x, Fe, t1)
stop = time()
print(f'Temps de calcul func basi : {stop-start}')
sf.write(f'results\\echo{t1}.flac', y, Fe)

start = time()
y = echo_scipy(x, Fe, t1)
stop = time()
print(f'Temps de calcul func scipy : {stop-start}')  

from math import pi
import scipy

n=np.arange(len(y)) #le temps discret que dure le signal y
f0=1261 #fréquence en Hertz
#La fréquence réduite est f0/Fe = 0.0263
parasite = 0.1*np.cos(2*pi*f0/Fe*n)
yb=y + parasite #on ajoute une onde parasite
sf.write(f'results\\echo_pollution{f0}.flac', yb, Fe)

def rejette1(entree,f0,Fs):

#valeur recommandée pour rho: au moins 0.9 et jamais plus de 1
    lfilter=scipy.signal.lfilter
    z0=np.exp(2*1j*pi*f0/Fs) # position sur le cercle unite de la fréquence à éliminer
    
    sortie=lfilter([1 ,-z0],1,entree)
    
    sortie=lfilter([1, -np.conj(z0)],1,sortie) #Le signal est reel, il faut donc éliminer -f0 aussi!
    
    
    sortie=np.real(sortie);# En raison d'erreurs d'arrondi il peut subsister une partie imaginaire
    return sortie
def rejette2(entree,f0,Fs,rho):

#valeur recommandée pour rho: au moins 0.9 et jamais plus de 1
    lfilter=scipy.signal.lfilter
    z0=np.exp(2*1j*pi*f0/Fs) # position sur le cercle unite de la fréquence à éliminer
    
    sortie=lfilter([1 ,-z0],[1 ,-rho*z0],entree)
    
    sortie=lfilter([1, -np.conj(z0)],[1 ,-rho*np.conj(z0)],sortie) #Le signal est reel, il faut donc éliminer -f0 aussi!

    sortie=np.real(sortie)# En raison d'erreurs d'arrondi il peut subsister une partie imaginaire
    return sortie

premier_filtre = rejette1(yb, f0, Fe)
second_filtre = rejette2(yb, f0+6, Fe, 0.90)

sf.write(f'results\\rejette1_{f0}.flac', premier_filtre, Fe) # rejette1 atténue beaucoup trop le signal que l'on doit garder
sf.write(f'results\\rejette2_{f0}.flac', second_filtre, Fe)

#rho < 1 pour la stabilité  

# Prendre un rho très proche de 1 va faire un filtre plus sélectif, mais il faut avoir une information
# très précise sur la fréquence à rejeter

