M=1:180;
r=0.83;
q=0.1;
P=25;

% 1.1
N=100000; % minták száma
u=zeros(N,1); % helyfoglalás a gerjeszt?jelnek
phi=rand(1,max(M)); % véletlen fázis generálás
for n=0:N-1
u(n+1)=sum(sin(2*pi*(M*n/N+phi)));
end

figure
plot(0:N-1,u);
title('u(n)');

% 1.2
[A,B]=analyse(u,max(M));
figure
subplot(2,1,1);
stem(A); % u(n) amplitúdóspektruma
title('u(n) spektruma az M. harmonikus frekvenciájáig');
subplot(2,1,2);
stem(B); % u(n) fázisspektruma


% 1.3
Hz=tf([(1-r) 0],[1 0 r],1); % átviteli függvény képzése
y=lsim(Hz,u); % az átviteli függvény gerjesztése u(n)-el

[Ay,By]=analyse(y,max(M)); % a kimen? jel amplitúdó- és fázisspektruma
H=Ay./A; % er?sítés
dP=By-B; % fázistolás

figure
subplot(2,1,1);
stem(H); 
title('Hz mért átvitele');
subplot(2,1,2);
stem(dP); 


% 1.4
% a 100. mintától már b?ven állandósult állapotban van y(n)
K=100;

X=toeplitz(zeros(P,1),[0 u(1:end-1)']); % regressziós vektor el?állítása
% a gerjeszt? jel elejére bekerült egy 0, így egy mintával késleltetve lett
% tehát nem kell figyelni kés?bb az y és u indexelésére
X=X';
x=X(K:end,:); % (K-1)-ig tartó mintasorozat eldobása

w0=ones(P,1)/P; % a súlyok inicializálása
P0=eye(P); % a P mátrix iniciílizálsa
w=recLS(w0,P0,x,y(K:end)); % rekurrzív LS függvény, w: súlytényez?k
y2=x*w;
figure
plot(0:9999,y2(1:10000),'b',0:9999,y(K:K+9999),'r');
legend('lin.kombinátor','adaptálandó rendszer');
title('Az adaptálandó rendszer és a lin. kombinátor kiemenete');

% 1.5
sor = (1-r)*(-r).^(1:P); % sorfejtett alak együtthatói
suly = impulse(Hz,1:2*P); % Hz impulzusválaszának els? 2P együtthatója 

figure
stem(sor);
hold on
stem(suly);
title('A sorfejtett alak együtthatói és Hz súlyfüggvénye');

figure
stem(w);
hold on
stem(suly);
title('A lineáris kombinátor együtthatói és Hz súlyfüggvénye');

y3=X*w; % lin. kombiátor kimenete a teljes u(n) gerjesztésre
[Alin,Blin]=analyse(y3,max(M)); % a lin.kombinátor kimenetének spektruma
figure
subplot(2,1,1);
% er?sítések
stem(Alin./A);
title('Az adaptálandó redszer és a lin. kombinátor er?sítése');
hold on
stem(Ay./A);
legend('lin.komb.','H(z)');
subplot(2,1,2);
stem(abs((Ay-Alin)./A)); % er?sítések különbsége


% fázistolások
figure
subplot(2,1,1);
stem(Blin-B);
title('Az adaptálandó redszer és a lin. kombinátor fázistolása');
hold on
stem(By-B);
legend('lin.komb.','H(z)');
subplot(2,1,2);
stem(abs(By-Blin)); % fázistolások különbsége

% 2 - LMS algoritmus

% a Wiener-Hopf egyenlet R mátrixának számítás - a bátorsági tényez?höz

R = X'*X;
R=R/max(M);
lam = max(eig(R)); % lambda az R mátrix sajátértékeinek maximuma
mu = 1/(10*lam); % bátorsági tényezõ

W = zeros(P,N); % súlymátrix inicializálása nullákkal

for i=P:N/2
    X=u(i-P+1:i);
    e=y(i)-X'*W(:,i);
    W(:,i+1)=W(:,i)+2*mu*X*e; % a súlyok rekurzív számítása LMS módszerrel
end

Hz_q = tf ([(1-(r-q)), 0, 0], [1, 0, 0, (r-q)], 1); % H(z) csökkentett r-el
y_q = lsim (Hz_q,u); % az új modell kimenete u(n) gerjesztésre

% a súlyok számításának folytatása az el?bbi módon
for i=N/2+1:N
    X=u(i-P+1:i);
    e=y_q(i)-X'*W(:,i);
    W(:,i+1)=W(:,i)+2*mu*X*e;
end

figure
% súlyok csökken? sorrendben
[maxW, maxNum] = sort(abs(W(:,end-1)),'descend');
hold all
for i=1:5
plot(abs(W(maxNum(i),:)')); % az 5 legnagyobb súly kirajzolása
end
legend(strcat('W',num2str(maxNum(1:5),-1)));
hold off
title('Az 5 legnagyobb súly konvergecia diagramja');

% 3

P2=6; % 6 együtthatót kell számolni
X=toeplitz(zeros(P2,1),[0 u(1:end-1)']);
R=X*X'/max(M); % ugyanúgy kiszámoljuk R-t a bátorsági tényez?höz
lam=max(eig(R));
mu=1/(10*lam); % bátorsági tényez?

% W1 és W2 számolása közös ciklusban
W1 = zeros(P2,N); % súlyok az normál modellhez
W2 = zeros(P2,N); % súlyok a csökkentett r-? modellhez

% equation-error formulation és LMS egyenletek
for i=P2:N-1
    X1=[u(i-1:i)', y(i-3:i)']; % 2 bemeneti és 4  korábbi kimeneti minta
    X2=[u(i-1:i)', y_q(i-3:i)'];
    e1=y(i+1)-X1*W1(:,i);
    e2=y_q(i+1)-X2*W2(:,i);
    W1(:,i+1)=W1(:,i)+2*mu*X1'*e1;
    W2(:,i+1)=W2(:,i)+2*mu*X2'*e2;
end

figure
for i=1:P2
    plot(W1(i,:)')
    hold on
end
legend('a1','a2','b1','b2','b3','b4');
title('Súlyok konvergencia diagramja a normál IIR rendszerhez');

figure

for i=1:P2
    plot(W2(i,:)')
    hold on
end
legend('a1','a2','b1','b2','b3','b4');
title('Súlyok konvergencia diagramja a csökkentett r-? IIR rendszerhez');



