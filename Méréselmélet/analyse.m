function [A,B]=analyse(a,M)
N=size(a,1); % minták száma
F=fft(a); % diszkrét fourier transzformáció
A=abs(F)*2/N; % amplitúdó spektrum
B=angle(F); % fázisspektrum
A=A(2:M+1); % az els? M frekvenciakomponens megtartása
B=B(2:M+1); 
end