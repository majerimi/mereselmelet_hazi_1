function [abc_ms,cov_ms] = MS(A,B,C,w)
% valtozok
global ma mb mc sa sb sc sw t0 dt f0 N;
% becslo
abc_ms=zeros(3,3); 
% variancia
cov_ms=zeros(3,3); 
% Sigma_aa mátrix
S_abc=diag([sa^2 sb^2 sc^2]); 
% a priori ismeret
m=[ma;mb;mc];

for i=1:3
   t=t0:dt:(t0+dt*(N(i)-1));
   % a veletlen parameterekbol osszeallitott jel
   U_t=(A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i)))';
   % megfigyelesi matrix
   U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(i))']; 
   % szumma_nn mátrix
   Sz_nn=diag(sw^2*ones(1,N(i))); 
   % explicit forma az a posteriori varhato ertekre
   abc_ms(i,:)=(m+(U'*Sz_nn^-1*U + S_abc^-1)^-1 * U'*Sz_nn^-1 * (U_t-U*m))';
   cov_ms(i,:)=diag((U'*Sz_nn^-1*U + S_abc^-1)^-1);
end
end