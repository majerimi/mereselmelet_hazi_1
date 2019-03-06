function Var = noiseVar(A,B,C,w,a)
t0=5;dt=0.9;f0=50/1000;
N=[3 30 300];
Var=zeros(3,1); % helyfoglalás a varianciáknak
for i=1:3
   n=0:1:(N(i)-1);
   t=t0+n*dt;
   z=A*sin(2*pi*f0*t)+B*cos(2*pi*f0*t)+C+w(1:N(i)); % jel generálása
   % megfigyelési mátrix
   U=a(i,1)*sin(2*pi*f0*t)+ a(i,2)*cos(2*pi*f0*t)+ a(i,3)*ones(1,N(i));
   Var(i) = var(z-U); % variancia számítása
end
end
