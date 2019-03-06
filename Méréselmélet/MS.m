function [p,v] = MS(A,B,C,w)
N=[3 10 100];
ma=1; mb=2; mc=0; sa=.1; sb=.2;
sc=.1; sw=.2;
t0=10;dt=3;f0=50/1000;
p=zeros(3,3); % helyfoglalás a becsl?knek
v=zeros(3,3); % helyfoglalás a varianciáknak
Sa=diag([sa^2 sb^2 sc^2]); % Sigma_aa mátrix
m=[ma;mb;mc];

for i=1:3
   n=0:1:(N(i)-1);
   t=t0+n*dt;
   % a véletéen paraméterekb?l összeállított jel
   z=(A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i)))';
   U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(i))']; % megfigyelési mátrix
   Sn=diag(sw^2*ones(1,N(i))); % Sigma_nn mátrix
   p(i,:)=(m+(U'*Sn^-1*U + Sa^-1)^-1 * U'*Sn^-1 * (z-U*m))';
   v(i,:)=diag((U'*Sn^-1*U + Sa^-1)^-1);
end
end