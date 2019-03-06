function [abc_ms,cov_ms] = ML(A,B,C,w)
% valtozok
global sw t0 dt f0 N;
% becslo
abc_ms=zeros(3,3);
% variancia
cov_ms=abc_ms;

for i=1:3
   t=t0:dt:(t0+dt*(N(i)-1));
   % a veletlen parameterekbol osszeallitott jel
   U_t=(A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i)))'; 
   % megfigyelesi matrix
   U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(i))'];
   Sz_nn=diag(sw^2*ones(1,N(i))); 
   % explicit forma az a posteriori varhato ertekre
   abc_ms(i,:)=((U'*Sz_nn^-1*U)^-1 * U'*Sz_nn^-1 *U_t)';
   cov_ms(i,:)=diag((U'*Sz_nn^-1*U)^-1)';
end
end 