function s = Synthesis(N,delta,f,a,phi,RSB);
% s = synthese(N,delta,f,a,phi);
% Synthèse d'un signal s de longueur N, constitué d'une somme de sinusoïdes amorties.

t = 0:N-1;
delta = min([delta(:),zeros(length(delta),1)],[],2);
logz = delta + i*2*pi*f(:);
alpha = a(:).*exp(i*phi(:));
x = sum( (alpha*ones(1,N)) .* exp(logz*t) ).';
if nargin<6,
	s = x;
else,
	Ex = real(x'*x) / N;
	b = randn(N,1) + i*randn(N,1);
	Eb = real(b'*b) / N;
	b = b * sqrt(Ex/Eb) * 10^(-RSB/20);
	s = x + b;
end;