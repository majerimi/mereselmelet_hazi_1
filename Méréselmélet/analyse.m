function [A,B]=analyse(a,M)
N=size(a,1); % mint�k sz�ma
F=fft(a); % diszkr�t fourier transzform�ci�
A=abs(F)*2/N; % amplit�d� spektrum
B=angle(F); % f�zisspektrum
A=A(2:M+1); % az els? M frekvenciakomponens megtart�sa
B=B(2:M+1); 
end