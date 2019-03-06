function [p,v] = ML(A,B,C,w)
sw=.1;
t0=5;dt=0.9;f0=50/1000;
N=[3 30 300];
% helyfoglal�s a becsl?knek �s varianci�knak
p=zeros(3,3);
v=p;

for i=1:3
   n=0:1:(N(i)-1);
   t=t0+n*dt;
   z=(A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i)))'; % a jel gener�l�sa
   U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(i))']; % megfigyel�si m�trix
   Sn=diag(sw^2*ones(1,N(i))); % Sigma_nn m�trix
   p(i,:)=((U'*Sn^-1*U)^-1 * U'*Sn^-1 *z)';
   v(i,:)=diag((U'*Sn^-1*U)^-1)';
end
