% Analyse d'un signal à l'aide de la TFCT
% Transformée de Fourier à Court Terme

clear all; close all;

% Paramètre d'interface   - Affichage du spectrogramme
affich = 1 ; % affich = 0 - Sans affichage
             % affich = 1 - Avec affichage
             % Note: cf. spectrogram sous matlab
             
% Paramètres du Vocodeur
M  = 1024;       % Longueur de la fenêtre
R  = M/4;        % Incrément des temps d'analyse, appelé hop size, t_a=pR
N  = 1024;       % Ordre de la tfd
w  = hanning(M); % Définition de la fenetre d'analyse de taille M
ws = w;          % Définition de la fenêtre de synthèse

if (N < M) 
    fprintf(2, 'N (%d) doit être plus grand que M (%d). \n\n', N, M); 
    return;
end

% Lecture du signal sonore et initialisation des paramètres associés
[x,Fe]  = audioread('salsa.wav');
x       = x(1:10*Fe, 1); x = x(:); % Signal audio, vecteur colonne monovoie
                                   % (voie gauche si stéréo)
Nx      = length(x);               % Longueur du signal en échantillons
Nt      = fix((Nx-M)/R);           % Nombre de trames

% Allocation mémoire
y = zeros(Nx, 1); % Signal de synthèse
if affich         % TFCT en cas d'affichage
    Xtilde = zeros(N, Nt);  
end

for p = 1:Nt              % Boucle sur les trames
   deb = (p-1)*R + 1;     % Début de trame
   fin = deb + M - 1;     % Fin de trame
   tx  = x(deb:fin) .* w; % Extraction de la trame pondérée par la fenêtre
   X   = fft(tx, N);      % Tfd à l'instant b
   
   if affich, Xtilde(:, p) = X; end
   
   % Début des opérations de transformation de TFCT X pour obtenir Y
   % .... (par défaut: pas de transformation, i.e. Y = X;
   
   Y = X;
   % Fin des opérations de transformation de TFCT X pour obtenir Y
   
   % Calcul de ys à partir de Y (à compléter)
   % .....
%   y(deb:fin) = y(deb:fin) + ys(1:M);   % (OLA) Addition recouvrement
end

%soundsc(y, Fe) % À décommenter pour jouer le son

% Affichage du spectrogramme
if affich    
    freq = Fe * (0:N/2)/N;            % Échelle des fréquences en Hz
    b    = ((M-1)/2 + R*(0:Nt-1))/Fe; % Instants d'analyse (au centre de la
                                      % fenêtre)
    imagesc(b, freq, db(abs(Xtilde(1:(N/2+1), :))));
    axis xy;
end