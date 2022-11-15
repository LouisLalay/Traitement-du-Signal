#%% SECTION de definition et d'imports a lancer au debut
import numpy as np
import soundfile as sf

import matplotlib.pyplot as plt
import scipy 
from scipy.signal import lfilter
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
randn=np.random.randn
def norm(X):
    return ((X**2).sum())**0.5
#%%
def une_colonne_spectrogramme(u,M,N,n):
    """ 
    Renvoie une colonne de spectrogramme c'est a dire la TFD de taille M
    d'une morceau de u debutant en n et de taille N multiplie par une fenêtre
     de hamming """
    uu=u.copy().reshape(-1)  # on recopie sinon on risque d'écraser 
    # construction de la fenêtre
    idx=np.arange(0,N)
    w=0.54-0.46*cos(2*pi*idx/(N-1))
    # les index tels que u(m)w(n-m) non nul
    m=np.arange(n,n+N).astype(np.int)
    morceauu=uu[m]  #morceau utile de u
    fenu=morceauu*w  
    Uc=fft(fenu,M) # on calcule la TFD 
    Uc=abs(Uc) # on s'intéresse seulement au module
    Uc=Uc[0:M//2+1] # pour un signal reel il suffit de garder la moitié
    return Uc 

def affiche_spectrogramme(u,N,M=None,nb=None,Fe=8192):
    """Affiche le specrogramme du signal u
     La taille des fenêtres est N
     Si M n'est pas fourni il est pris égal à N
     nb est le pas entre deux fenêtres dont on calcule la TFD 
     si nb n'est pas fourni, nb est pris égal a N/2"""
    
    if M is None:
        M=N
    if nb is None:
        nb=N/2
    # On commence par créer un tableau de la bonne taille don les colonnes seront
    # calculées par une_colonne_spectrogramme
    uu=u.copy().reshape(-1)
    L=len(u)
    nombre_fen=int((L-N)//nb+1)
    spectro=np.zeros((M//2+1,nombre_fen))
    for k in range(nombre_fen):
        spectro[:,k]=une_colonne_spectrogramme(u,M,N,k*nb)
    temps_debut=0
    temps_fin=1/Fe*nb*nombre_fen 
    freq_debut=0
    freq_fin=(M/2)*(1/M)*Fe
    
    # ci-dessous de la cuisine python pour un affichage d'une image deux fois 
    # pus large que haute
    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(np.log(np.flipud(spectro)), interpolation='none', # cmap=plt.cm.Reds
              extent=[temps_debut,temps_fin,freq_debut,freq_fin])
    ax.set_aspect(1/2*(temps_fin/freq_fin))
    return spectro

def moindres_carres(p,x,z):
    """ % Cette fonction renvoie simplement le filtre h qui minimise l'énergie de 
% z-h*x    (h*x est la convolution PAS LE PRODUIT de h et x)
% h est de taille p+1
% Si la taille de z n'est pas suffisante, le programme se permet de tronquer
% les signaux en conséquence.
% x et z sont traites comme des signaux causaux
""" 
    assert (len(x.shape)==2 and min(x.shape)==1) or len(x.shape)==1, "La taille de x est mauvaise"
    assert (len(z.shape)==2 and min(z.shape)==1) or len(z.shape)==1, "La taille de z est mauvaise"
    xc=x.copy().reshape(-1) # colonne
    zc=z.copy().reshape(-1)
    # """ xmat est une matrice telle que h@xmat est h*x (h convolué avec x)"""
    xmat=np.zeros((p+1,len(xc)+p))
    for k in range(p+1):
        xmat[k,k:k+len(x)]=xc
    
    # si z est trop long
    if len(zc)>len(xc)+p:
        zc=zc[:len(xc)+p]
    # si z est trop court
    if len(zc)<len(xc)+p:
        zc=np.concatenate((zc,np.zeros(len(xc)+p-len(zc))))
    # resolution du problème : trouver h telle que norm(z-h@xmat) soit minimale
    # Vous pouvez essayer de démontrer la formule hors du TP
    L=xmat@(xmat.T)
    V=xmat@zc.T
    h=np.linalg.inv(L)@V
    return h

#%% PARTIE 2
h=np.asarray([1,1/2])
N=1000
X=randn(N)
sigma=0.1
Y=sigma*randn(N)

Z=lfilter(h,[1],X)+Y
hestime=moindres_carres(1,X,Z)
hestime

        
#%%PARTIE 3
#lecture des données 
fich=scipy.io.loadmat('parasite.mat')
Xparasite=fich['Xparasite'].reshape(-1)
Zparasite=fich['Zparasite'].reshape(-1)
Fe=22050
play(Zparasite,Fe)
#%%

hestim=moindres_carres(p_inconnu,Xparasite,Zparasite) #COMPLETER p_inconnu
Y=Zparasite-lfilter(hestim,1,Xparasite)
play(Y,Fe)

#%% PARTIE 4
def lpc_morceau(p,x):
    """ % Renvoie un filtre de taille p+1 dont la première valeur est 1 et qui
% minimise le h*x
""" 
    assert (len(x.shape)==2 and min(x.shape)==1) or len(x.shape)==1, "La taille de x est mauvaise"
    xc=x.copy().reshape(-1) # colonne
    # """ xmat est une matrice telle que h@xmat est h*x (h convolué avec x)"""
    xmat=np.zeros((p+1,len(xc)+p))
    for k in range(p+1):
        xmat[k,k:k+len(x)]=xc
    #n=np.arange(0,len(xc))
    #xc=xc*cos(pi*n/len(xc))
    
    # resolution du problème : trouver h telle que norm(z-h@xmat) soit minimale
    # Vous pouvez essayer de démontrer la formule hors du TP
    M=xmat@(xmat.T)
    L=M[1:,1:]
    V=M[1:p+1,0]
    h1=-np.linalg.inv(L)@V
    h=np.concatenate((np.ones(1),h1))
    return h

#%%
N=1000
h=np.asarray([1,1/2])
sigma=1
epsilon=sigma*randn(N)
X=lfilter([1],h,epsilon)
hestim=lpc_morceau(1,X)
hestim
#%% PARTIE 4.2
def lpc(signal,p,nb):
    """ fonction qui renvoie les meilleurs coefficients de prediction linéaire 
       pour les morceaux de signal de taille nb. Il y a, en gros, N/nb lignes
       ou N est le nombre d'échantillons du signal.
       Renvoie également la puissance du résidu epsilon"""
      
    sc=signal.copy().reshape(-1)
    sc=lfilter([1,-0.95],1,sc) #un filtrage passe haut pour renforcer les hautes fréquences
    N=len(sc)
    r=N%nb
    if r==0:
        taille=N//nb
    else:
        taille=N//nb+1
    out=np.zeros((taille,p+1))
    res=np.zeros(taille)
    for k in range(taille):
        deb=k*nb
        fin=min((k+1)*nb,len(sc))
        tmp=sc[deb:fin]
        
        if len(tmp)<nb:
            tmp=np.concatenate((tmp,np.zeros(nb-len(tmp))))
        
        h=lpc_morceau(p,tmp)
        out[k,:]=h
        epsilon=lfilter(h,1,tmp)
        
        res[k]=((epsilon**2).sum()/nb)**0.5
    return out,res

def joue_lpc(lpcs,res,nb):
    """joue le resultat de lpcs: cree des trames de bruit, les filtre par les coefficients de la LPC et renvoie un signal concaténé."""
    taille=lpcs.shape[0]
    out=np.zeros(taille*nb)
    cordesvocales=randn(len(out))
    # Alternative les cordes vocales envoient des impulsions régulières
    #cordesvocales=np.zeros(len(out))
    #ordesvocales[::NOMBRE_inconnu]=1   #COMPLETER 
    for k in range(taille):
        epsilon=res[k]*cordesvocales[k*nb:(k+1)*nb]
        
        tmp=lfilter([1],lpcs[k,:],epsilon)
        out[k*nb:(k+1)*nb]=tmp

    return out
#%%          
Fe=22050
fich=scipy.io.loadmat('modelear.mat')
aeiou=fich['aeiou']

lpcs,res=lpc(aeiou,15,2000) # tester d'autres valeurs de 15
out=joue_lpc(lpcs,res,2000)
play(out,Fe)
#%%
joue_lpc(lpcs,res)

#%% BONUS: Elimination des petits clicks (il faut utiliser des impulsions synchrones 
# dans les deux calculs de out0 et out1, par exemple avec des impulsions tous les 200
# échantillonns        
def multiplie_triangles(sig,nb):
    taille=int(len(sig)/nb)
    triangle=(nb/2-abs(np.arange(nb)-nb/2))/nb*2
    out=sig.copy()
    for k in range(0,taille,nb):
        out[k:k+nb]=out[k:k+nb]*triangle
    return out

nb=2000
lpcs0,res0=lpc(aeiou,15,nb)
lpcs1,res1=lpc(np.concatenate((aeiou[nb//2:],np.zeros((nb//2,1))),axis=0),15,nb)
out0=joue_lpc(lpcs0,res0,nb)
out1=joue_lpc(lpcs1,res1,nb)
out0=multiplie_triangles(out0,nb)
out1=multiplie_triangles(out1,nb)
out=out0.copy()
out[nb//2:]=out[nb//2:]+out1[:-nb//2]
play(out,Fe)