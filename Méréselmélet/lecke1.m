clc
clear all
r=0.87;
q=0.13;
M=250;
P=19;

%%%%% 1. feladat %%%%%
for turn=1:2 % P=1*P illetve 2*P esetek
    Pi = turn*P; % aktu�lis P be�ll�t�sa
    db=0:1:M-1; % szinusz mintasz�m vektor
    db=db/M; % szinusz mintasz�m vektor norm�l�sa

    s=sin(2*pi*db); % alapharm�nikus elk�sz�t�se
    phase = 2*pi*rand(M/2-1,1); % v�ltelen kezd�f�zisuk gener�l�sa a felharm�nikusoknak

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1)); % felharm�nikusok gener�l�sa, majd hozz�ad�sa az alapharm�nikushoz
    end

    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1); % modellezend� rendszer �tviteli f�ggv�nye
    y = lsim (system, s); % sz�rt mint�k

    fprintf('dimW=%dP eset �tviteli f�ggv�ny egy�tthat�i',turn);
    filt_ser=impulse(system,(1:Pi)) % �tviteli f�ggv�ny egy�tthat�i

    top = [0 s(1:end-1)]; % X m�trix els� sora
    left = zeros(Pi,1); % X m�trix els� oszlopa
    X   = toeplitz(left,top); % X m�trix el��ll�t�sa

    R = X*X'; % R m�trix kisz�m�t�sa
    R=R/M; % R m�trix norm�l�sa

    p = X * y; % P m�trix kisz�m�t�sa
    p = p/M; % P m�trix norm�l�sa

    W = R\p; % line�ris kombin�tor s�lyt�nyez�i, inv(R)*P
       
    fprintf('dimW=%dP eset egy�tthat�i',turn);
    
    fprintf('dimW=%dP eset sorfejtett alak �s a line�ris kombin�tor egy�tthat�inak elt�r�se',turn);
    filt_ser(1:Pi)-W;
    
    sys = tf([0; W]',1,1,'variable','z^-1'); % modellezett rendszer el��ll�t�sa
    ykalap = lsim (sys, s); % becs�lt mint�k gener�l�sa
    
    %%%Diagrammok elk�sz�t�se
    
    figure(1+turn-1)
    subplot(2,1,1)
    stem((1:Pi), filt_ser(1:Pi))
    hold on
    stem((1:Pi),W,'r')
    hold off
    subplot(2,1,2)
    stem((1:Pi), filt_ser(1:Pi)-W)
    
    figure(3+turn-1)
    subplot(2,1,1)
    stem((-M/2:M/2-1),abs(fft(y)/M))
    hold on
    stem((-M/2:M/2-1),abs(fft(ykalap)/M),'r')
    hold off
    subplot(2,1,2)
    stem((-M/2:M/2-1),abs(fft(y)/M)-abs(fft(ykalap)/M))

    figure(5+turn-1)
    subplot(2,1,1)
    stem(wrapToPi(angle(fft(y)/M)))
    ylim([-pi pi])
    hold on
    stem(wrapToPi(angle(fft(ykalap)/M)),'r')
    ylim([-pi pi])
    hold off
    subplot(2,1,2)
    stem(wrapToPi(wrapToPi(angle(fft(y)/M))-wrapToPi(angle(fft(ykalap)/M))))
    ylim([-pi pi])
end

%%%%%% 2. feladat %%%%%%%

for turn=1:2
    Pi = turn*P;
    db=0:1:M-1;
    db=db/M;

    s=sin(2*pi*db);
    phase = 2*pi*rand(M/2-1,1);

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1));
    end

    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
    y = lsim (system, s);

    filt_ser=impulse(system,(1:Pi));

    top = [0 s(1:end-1)];
    left = zeros(Pi,1);

    X   = toeplitz(left,top);

    R = X*X';
    R=R/M;

    lambda = max(eig(R)); % lambda az R m�trix saj�t�rt�keinek maximuma
    LMS_mu = 1/(3000*lambda); % b�tors�gi t�nyez� defini�l�sa
    szinusz = s;

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1));
    end

    for i=1:99
        szinusz = [szinusz s]; % 25000 mint�b�l �ll� esem�nyt�r l�trehoz�sa
    end
    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
    y = lsim (system, szinusz);

    W = zeros(Pi,25000); % kezdeti m�trix null�z�sa

    for i=Pi:25000/2
        X=szinusz(i-Pi+1:i); % P darab minta kiv�tele
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*X'*e; % LMS algoritmus sz�mol�sa
    end

    system = tf ([(1-(r-q)), 0, 0], [1, 0, 0, (r-q)], 1); % r q-val val� cs�kkent�se
    y = lsim (system, szinusz);

    for i=25000/2+1:25000
        X=szinusz(i-Pi+1:i); % P darab minta kiv�tele
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*X'*e; % LMS algoritmus sz�mol�sa
    end

    fprintf('dimW=%dP eset egy�tthat�i, eredeti',turn);
    W(:,25000/2)
    
    fprintf('dimW=%dP eset egy�tthat�i, cs�kkentett r',turn);
    W(:,25000)
    
    figure(7+turn-1)
    [maxW, maxNum] = sort(abs(W(:,end-1)),'descend'); % erem�nyek cs�kken� sorrendbe rendez�se
    hold all
    for i=1:5
        plot(abs(W(maxNum(i),:)')); % 5 legnagyobb egy�tthat� kirajzol�sa
    end
    legend(strcat('W',num2str(maxNum(1:5),-1)))
    hold off
end

%%%%% 3. feladat %%%%%

for turn=1:2
    Pi = turn*P;
    db=0:1:M-1;
    db=db/M;

    s=sin(2*pi*db);
    phase = 2*pi*rand(M/2-1,1);

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1));
    end

    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
    y = lsim (system, s);

    filt_ser=impulse(system,(1:Pi));

    top = [0 s(1:end-1)];
    left = zeros(Pi,1);

    X   = toeplitz(left,top);

    R = X*X';
    R=R/M;

    lambda = max(eig(R));
    LMS_mu = 1/(10*lambda);
    szinusz = s;

    lambda = 0.9; % feladat �lltal megadott param�terek
    v = 0.1;

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1));
    end

    for i=1:99
        szinusz = [szinusz s]; % 25000 mint�b�l �ll� esem�nyt�r l�trehoz�sa
    end
    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
    y = lsim (system, szinusz);

    W = zeros(Pi,25000);
    R = eye(Pi); % kezdeti m�trixok null�z�sa

    for i=Pi:25000/2
        X=szinusz(i-Pi+1:i); % P darab minta kiv�tele
        num = (R*X'*X*R);
        den = ((lambda/v)+X*R*X');
        R = (1/lambda)*(R-(num/den));
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*R*X'*e; % algoritmus sz�mol�sa
    end

    system = tf ([(1-(r-q)), 0, 0], [1, 0, 0, (r-q)], 1); % r q-val val� cs�kkent�se
    y = lsim (system, szinusz);

    for i=25000/2+1:25000
        X=szinusz(i-Pi+1:i); % P darab minta kiv�tele
        num = (R*X'*X*R);
        den = ((lambda/v)+X*R*X');
        R = (1/lambda)*(R-(num/den));
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*R*X'*e; % algoritmus sz�mol�sa
    end
    
    fprintf('dimW=%dP eset egy�tthat�i, eredeti',turn);
    W(:,25000/2)
    
    fprintf('dimW=%dP eset egy�tthat�i, cs�kkentett r',turn);
    W(:,25000)
    
    figure(9+turn-1)
    [maxW, maxNum] = sort(abs(W(:,end-1)),'descend');
    hold all
    for i=1:5
        plot(abs(W(maxNum(i),:)')); % 5 legnagyobb egy�tthat� kirajzol�sa
    end
    legend(strcat('W',num2str(maxNum(1:5),-1)))
    hold off
end

%%%%%% 4. feladat %%%%%%%
P=6; % 6 egy�tthat�t kell keresni

db=0:1:M-1;
db=db/M;

s=sin(2*pi*db);
phase = 2*pi*rand(M/2-1,1);

for i=2:M/2
    s = s+sin(i*2*pi*db+phase(i-1));
end

system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
y = lsim (system, s);

filt_ser=impulse(system,(1:P));

top = [0 s(1:end-1)];
left = zeros(P,1);

X = toeplitz(left,top);

R = X*X';
R=R/M;

lambda = max(eig(R));
LMS_mu = 1/(1000*lambda); % b�tors�gi t�nyez� meghat�roz�sa
szinusz = s;

for i=2:M/2
    s = s+sin(i*2*pi*db+phase(i-1));
end

for i=1:99
    szinusz = [szinusz s];
end
system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
y = lsim (system, szinusz);

W = zeros(P,25000); % kezdeti m�trix null�z�sa

for i=P:25000-1
    X=[szinusz(i-1:i), y(i-3:i)']; % a megadott formula szerint 2 gerjeszt� minta,
                                   % valamint 4 kimeneti minta kiv�laszt�sa
    e=y(i+1)-X*W(:,i);
    W(:,i+1)=W(:,i)+2*LMS_mu*X'*e;
end

fprintf('Egy�tthat�i:',turn);
W(:,end)

figure(11)
hold all
for i=1:6
	plot(W(i,:)');
end
legend('W0','W1','W2','W3','W4','W5',-1)
plot(W')
