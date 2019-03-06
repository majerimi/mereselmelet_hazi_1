function [variancie,varvariancie] = noise_variancie(A,B,C,w,a)
global t0 dt f0 N;

variancie=zeros(3,1);
varvariancie=zeros(3,1);

for i=1:3
   t=t0:dt:(t0+dt*(N(i)-1));
   % a veletlen parameterekbol osszeallitott jel
   U_t=A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i)); 
   % megfigyelesi matrix
   U=a(i,1)*sin(2*pi*f0*t)+ a(i,2)*cos(2*pi*f0*t)+ a(i,3)*ones(1,N(i));
   % variancia szamitas
   variancie(i) = var(U_t-U);
   varvariancie(i) = var(var(U_t-U));
end
end
