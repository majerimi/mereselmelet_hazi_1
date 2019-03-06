function p = LS(A,B,C,w)
t0=5;dt=0.9;f0=50/1000;
N=[3 30 300];
p=zeros(3,3);


for i=1:3
   n=0:1:(N(i)-1);
   t=t0+n*dt;
   z=A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i));
   U=[sin(2*pi*f0*t)' cos(2*pi*f0*t)' ones(1,N(i))'];
   p(i,:)=(U'*U)^-1 * U' * z';
end
end