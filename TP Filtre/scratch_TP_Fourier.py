#%% SECTION de definition et d'imports a lancer au debut
import numpy as np
import soundfile as sf

import matplotlib.pyplot as plt
import scipy 
import scipy.io
import platform
import time

import tempfile
import os

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
        print ('Vous allez cr�er un fichier son de plus de 20 secondes.')
        rep=None
        while rep is None:
            x=input('Voulez-vous continuer? (o/n)')
            if x=='o':
                rep=True
            if x=='n':
                rep=False
            if rep is None:
                print ('R�pondre par o ou n, merci. ')
    if rep:
        sf.write(fichier,z,Fe)
        os.system('/usr/bin/play '+fichier+' &')



#%% quelques simplifications de fonctions usuelles
exp=np.exp
cos=np.cos

fft=np.fft.fft
ifft=np.fft.ifft
real=np.real

plot=plt.plot
stem=plt.stem
show=plt.show # force l'affichage du graphique courant
i=np.complex(0,1)
pi=np.pi 

#%% UNE ONDE ET SA FFT

nu=0.123
n=np.arange(0,100)
onde=np.exp(2*i*pi*nu*n) # est-ce vraiment une onde de Fourier de la TFD?
stem(n,abs(fft(onde)))
show()
stem(n,abs(fft(real(onde))))

#%%
a=np.asarray([1,2,3])
b=np.asarray([1,1,0])
c=np.asarray([4,3,5]) #COMPLETER convolution des deux a calculer a la main (convolution circulaire taille 3)

fft(a)*fft(b)-fft(c)
#%%
ifft(fft([1,1,1]))
ifft(fft([1,1,1,0]))
ifft(fft([1,2,3])*fft([1,1,0]))-[4,3,5]
#%%
h=np.asarray([1,2,3])
g=np.asarray([1,2,1])
N=max(len(h),len(g))
hz=np.concatenate((h,np.zeros((N-len(h))))) #COMPLETER
gz=np.concatenate((g, np.zeros((N-len(g))))) #COMPLETER
out=real(ifft(fft(hz)*fft(gz)))
out=out[0:len(h)+len(g)-1]

def convolution_rapide(g,h):
    N=max(len(g),len(h)) 
    hz=np.concatenate((h, np.zeros((N-len(h))))) #COMPLETER
    gz=np.concatenate((g, np.zeros((N-len(g))))) #COMPLETER
    result=real(ifft(fft(hz)*fft(gz)))
    return result

h=np.random.randn(10000)
g=np.random.randn(10000)
t0=time.time()         
res=convolution_rapide(h,g)
print ('temps rapide=',time.time()-t0)

t0=time.time()         
res=scipy.convolve(h,g)
print ('temps normal=',time.time()-t0)

#%%
Fe=1000 ;Te=1/Fe ; #�chantillonnage a 1000 Hz 
t=np.arange(0,300)*Te# t repr�sente le temps en secondes de 300 �chantillons

xi=115; # fr�quence de l'onde 
h=cos(2*pi*xi*t) # equiv. cos(2*pi*(xi/Fe)*np.arange(0,300))
fh=fft(h)
fh=fh[0:151]
plot(np.arange(0,151)/300*Fe,abs(fh))

#%% une onde
exp=np.exp
i=np.complex(0,1)
pi=np.pi 
n=np.arange(0,10000)
f=0.01
u=exp(2*i*pi*f * n)
play(np.real(u),44100)

