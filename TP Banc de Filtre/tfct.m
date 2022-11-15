% analyse d'un signal à l'aide de la TFCT
% transformée de Fourier à Court Terme
clear all; close all;

[x,Fe] = wavread('song.wav');
x = x(1:10*Fe,1);x = x(:); % ainsi x est un vect. colonne
        % monovoie (voie gauche si stéréo)

N = length(x); % longueur du signal
Nw = 512;
w = hanning(Nw); % définition de la fenetre d'analyse
ws = w; % définition de la fenêtre de synthèse
R = Nw/4; % incrément sur les temps d'analyse,
          % appelé hop size, t_a=uR
M = 512; % ordre de la tfd
L = M/2+1;

affich = 1 ; % pour affichage du spectrogramme, 0 pour
             % pour faire analyse/modif/synthèse sans affichage
             % note: cf. spectrogram sous matlab


Nt = fix( (N-Nw)/R ); % calcul du nombre de tfd à calculer
y = zeros(N,1); % signal de synthèse

if affich
    Xtilde = zeros(M,Nt);
end

for u=1:Nt;  % boucle sur les trames
   deb = (u-1)*R +1; % début de trame
   fin = deb + Nw -1; % fin de trame
   tx = x(deb:fin).*w; % calcul de la trame  
   X = fft(tx,M); % tfd à l'instant b
   
   if affich, Xtilde(:,u)=X;end
   
   % opérations de transformation (sur la partie \nu > 0)
   % ....
   Y = X;
   % fin des opération de transformation
   
   % resynthèse
   % overlap add
end



soundsc(y,Fe)

if affich
    
freq = (0:M/2)/M*Fe;
b = [0:Nt-1]*R/Fe+Nw/2;
imagesc(b,freq,db(abs(Xtilde(1:L,:))))
axis xy

end