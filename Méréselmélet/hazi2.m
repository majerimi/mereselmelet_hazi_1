M=1:180;
r=0.83;
q=0.1;
P=25;

% 1.1
N=100000; % mint�k sz�ma
u=zeros(N,1); % helyfoglal�s a gerjeszt?jelnek
phi=rand(1,max(M)); % v�letlen f�zis gener�l�s
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
stem(A); % u(n) amplit�d�spektruma
title('u(n) spektruma az M. harmonikus frekvenci�j�ig');
subplot(2,1,2);
stem(B); % u(n) f�zisspektruma


% 1.3
Hz=tf([(1-r) 0],[1 0 r],1); % �tviteli f�ggv�ny k�pz�se
y=lsim(Hz,u); % az �tviteli f�ggv�ny gerjeszt�se u(n)-el

[Ay,By]=analyse(y,max(M)); % a kimen? jel amplit�d�- �s f�zisspektruma
H=Ay./A; % er?s�t�s
dP=By-B; % f�zistol�s

figure
subplot(2,1,1);
stem(H); 
title('Hz m�rt �tvitele');
subplot(2,1,2);
stem(dP); 


% 1.4
% a 100. mint�t�l m�r b?ven �lland�sult �llapotban van y(n)
K=100;

X=toeplitz(zeros(P,1),[0 u(1:end-1)']); % regresszi�s vektor el?�ll�t�sa
% a gerjeszt? jel elej�re beker�lt egy 0, �gy egy mint�val k�sleltetve lett
% teh�t nem kell figyelni k�s?bb az y �s u indexel�s�re
X=X';
x=X(K:end,:); % (K-1)-ig tart� mintasorozat eldob�sa

w0=ones(P,1)/P; % a s�lyok inicializ�l�sa
P0=eye(P); % a P m�trix inici�liz�lsa
w=recLS(w0,P0,x,y(K:end)); % rekurrz�v LS f�ggv�ny, w: s�lyt�nyez?k
y2=x*w;
figure
plot(0:9999,y2(1:10000),'b',0:9999,y(K:K+9999),'r');
legend('lin.kombin�tor','adapt�land� rendszer');
title('Az adapt�land� rendszer �s a lin. kombin�tor kiemenete');

% 1.5
sor = (1-r)*(-r).^(1:P); % sorfejtett alak egy�tthat�i
suly = impulse(Hz,1:2*P); % Hz impulzusv�lasz�nak els? 2P egy�tthat�ja 

figure
stem(sor);
hold on
stem(suly);
title('A sorfejtett alak egy�tthat�i �s Hz s�lyf�ggv�nye');

figure
stem(w);
hold on
stem(suly);
title('A line�ris kombin�tor egy�tthat�i �s Hz s�lyf�ggv�nye');

y3=X*w; % lin. kombi�tor kimenete a teljes u(n) gerjeszt�sre
[Alin,Blin]=analyse(y3,max(M)); % a lin.kombin�tor kimenet�nek spektruma
figure
subplot(2,1,1);
% er?s�t�sek
stem(Alin./A);
title('Az adapt�land� redszer �s a lin. kombin�tor er?s�t�se');
hold on
stem(Ay./A);
legend('lin.komb.','H(z)');
subplot(2,1,2);
stem(abs((Ay-Alin)./A)); % er?s�t�sek k�l�nbs�ge


% f�zistol�sok
figure
subplot(2,1,1);
stem(Blin-B);
title('Az adapt�land� redszer �s a lin. kombin�tor f�zistol�sa');
hold on
stem(By-B);
legend('lin.komb.','H(z)');
subplot(2,1,2);
stem(abs(By-Blin)); % f�zistol�sok k�l�nbs�ge

% 2 - LMS algoritmus

% a Wiener-Hopf egyenlet R m�trix�nak sz�m�t�s - a b�tors�gi t�nyez?h�z

R = X'*X;
R=R/max(M);
lam = max(eig(R)); % lambda az R m�trix saj�t�rt�keinek maximuma
mu = 1/(10*lam); % b�tors�gi t�nyez�

W = zeros(P,N); % s�lym�trix inicializ�l�sa null�kkal

for i=P:N/2
    X=u(i-P+1:i);
    e=y(i)-X'*W(:,i);
    W(:,i+1)=W(:,i)+2*mu*X*e; % a s�lyok rekurz�v sz�m�t�sa LMS m�dszerrel
end

Hz_q = tf ([(1-(r-q)), 0, 0], [1, 0, 0, (r-q)], 1); % H(z) cs�kkentett r-el
y_q = lsim (Hz_q,u); % az �j modell kimenete u(n) gerjeszt�sre

% a s�lyok sz�m�t�s�nak folytat�sa az el?bbi m�don
for i=N/2+1:N
    X=u(i-P+1:i);
    e=y_q(i)-X'*W(:,i);
    W(:,i+1)=W(:,i)+2*mu*X*e;
end

figure
% s�lyok cs�kken? sorrendben
[maxW, maxNum] = sort(abs(W(:,end-1)),'descend');
hold all
for i=1:5
plot(abs(W(maxNum(i),:)')); % az 5 legnagyobb s�ly kirajzol�sa
end
legend(strcat('W',num2str(maxNum(1:5),-1)));
hold off
title('Az 5 legnagyobb s�ly konvergecia diagramja');

% 3

P2=6; % 6 egy�tthat�t kell sz�molni
X=toeplitz(zeros(P2,1),[0 u(1:end-1)']);
R=X*X'/max(M); % ugyan�gy kisz�moljuk R-t a b�tors�gi t�nyez?h�z
lam=max(eig(R));
mu=1/(10*lam); % b�tors�gi t�nyez?

% W1 �s W2 sz�mol�sa k�z�s ciklusban
W1 = zeros(P2,N); % s�lyok az norm�l modellhez
W2 = zeros(P2,N); % s�lyok a cs�kkentett r-? modellhez

% equation-error formulation �s LMS egyenletek
for i=P2:N-1
    X1=[u(i-1:i)', y(i-3:i)']; % 2 bemeneti �s 4  kor�bbi kimeneti minta
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
title('S�lyok konvergencia diagramja a norm�l IIR rendszerhez');

figure

for i=1:P2
    plot(W2(i,:)')
    hold on
end
legend('a1','a2','b1','b2','b3','b4');
title('S�lyok konvergencia diagramja a cs�kkentett r-? IIR rendszerhez');



