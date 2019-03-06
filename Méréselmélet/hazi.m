M=100;
N=[3 30 300];
ma=2; mb=1; mc=0; sa=.2; sb=.1;
sc=.3; sw=.1;
t0=5;dt=0.9;f0=50/1000;
% Helyfoglal�s az 1.2 �s 1.5-�s feladatok v�ltoz�inak
A_ms=zeros(3,3,M);
A_ls=A_ms;

% V�letlen param�tereg gener�l�sa
A1=normrnd(ma,sa);
B1=normrnd(mb,sb);
C1=normrnd(mc,sc);
w1=normrnd(0,sw,[1,max(N)]);
[a_ms,var_ms] = MS(A1,B1,C1,w1); % 1.1
[a_ml,var_ml] = ML(A1,B1,C1,w1); % 1.3

A_ms(:,:,1)=a_ms;
A_ls(:,:,1)=a_ml; % az LS becsl? (1.4) megegyezik az ML becsl?vel

% 1.1 �s 1.4 m�r�sek megism�tl�se (M-1)-szer
for  k=2:M
    A=normrnd(ma,sa);
    B=normrnd(mb,sb);
    C=normrnd(mc,sc);
    w=normrnd(0,sw,[1,max(N)]);
    A_ms(:,:,k)=MS(A,B,C,w); % 1.2
    A_ls(:,:,k)=ML(A,B,C,w); % 1.5
end

% 1.2, MS becslok �tlaga �s varianci�ja
meanA_ms = mean(A_ms,3);
varA_ms = var(A_ms,[],3);

% 1.5, LS becslok �tlaga �s varianci�ja
meanA_ls = mean(A_ls,3);
varA_ls = var(A_ls,[],3);

% 1.7, a csatornazaj varianci�ja
var_noise = noiseVar(A1,B1,C1,w1, a_ms); 

% 1.8, rekurz�v LS becsl? sz�m�t�sa
a0=A_ls(1,:)';
n=0:1:(N(3)-1);
t=t0+n*dt;
U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(3))'];
z=(A1*sin(2*pi*f0*t)+B1*cos(2*pi*f0*t)+C1+w1)';
P0=(U(1:N(1))'*U(1:N(1)))^-1;
[a_recLS, G, e] = recLS(a0,P0,U,z);

% v�ltoz�k megjelen�t�se
figure
plot(n(1:end-1)*dt,G(1,:),'r',n(1:end-1)*dt,G(2,:),'g',n(1:end-1)*dt,G(3,:),'b')
grid on
figure
plot(n(1:end-1)*dt,G(1,:),'r',n(1:end-1)*dt,G(2,:),'g',n(1:end-1)*dt,G(3,:),'b')
grid on
figure
plot(n(1:end-1)*dt,e,'r')
grid on
