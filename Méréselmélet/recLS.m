function [a ,G, e] = recLS(a0,P0,u,z)
N=900;

% kezdeti értékek beállítása
a=a0;
P=P0;
% rekurzív LS becsl? egyenletei
for i=1:N
   G=P*u(i,:)'/(1+u(i,:)*P*u(i,:)');
   e=z(i)-u(i,:)*a;
   a=a+G*e;
   P=(eye(25)-G*u(i,:))*P;
end
end