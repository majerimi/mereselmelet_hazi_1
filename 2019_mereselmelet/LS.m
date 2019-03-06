function abc_ls = LS(A,B,C,w)
global t0 dt f0 N;

abc_ls=zeros(3,3);


for i=1:3
   t=t0:dt:(t0+dt*(N(i)-1));
   z=A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i));
   U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(i))'];
   abc_ls(i,:)=(U'*U)^-1 * U' * z';
end
end