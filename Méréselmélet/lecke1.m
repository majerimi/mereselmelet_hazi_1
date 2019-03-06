clc
clear all
r=0.87;
q=0.13;
M=250;
P=19;

%%%%% 1. feladat %%%%%
for turn=1:2 % P=1*P illetve 2*P esetek
    Pi = turn*P; % aktuális P beállítása
    db=0:1:M-1; % szinusz mintaszám vektor
    db=db/M; % szinusz mintaszám vektor normálása

    s=sin(2*pi*db); % alapharmónikus elkészítése
    phase = 2*pi*rand(M/2-1,1); % véltelen kezdõfázisuk generálása a felharmónikusoknak

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1)); % felharmónikusok generálása, majd hozzáadása az alapharmónikushoz
    end

    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1); % modellezendõ rendszer átviteli függvénye
    y = lsim (system, s); % szûrt minták

    fprintf('dimW=%dP eset átviteli függvény együtthatói',turn);
    filt_ser=impulse(system,(1:Pi)) % átviteli függvény együtthatõi

    top = [0 s(1:end-1)]; % X mátrix elsõ sora
    left = zeros(Pi,1); % X mátrix elsõ oszlopa
    X   = toeplitz(left,top); % X mátrix elõûllítása

    R = X*X'; % R mátrix kiszámítása
    R=R/M; % R mátrix normálása

    p = X * y; % P mátrix kiszámítása
    p = p/M; % P mátrix normálása

    W = R\p; % lineáris kombinátor súlytényezõi, inv(R)*P
       
    fprintf('dimW=%dP eset együtthatói',turn);
    
    fprintf('dimW=%dP eset sorfejtett alak és a lineáris kombinátor együtthatóinak eltérése',turn);
    filt_ser(1:Pi)-W;
    
    sys = tf([0; W]',1,1,'variable','z^-1'); % modellezett rendszer elõállítása
    ykalap = lsim (sys, s); % becsült minták generálása
    
    %%%Diagrammok elkészítése
    
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

    lambda = max(eig(R)); % lambda az R mátrix sajátértékeinek maximuma
    LMS_mu = 1/(3000*lambda); % bátorsági tényezõ definiálása
    szinusz = s;

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1));
    end

    for i=1:99
        szinusz = [szinusz s]; % 25000 mintából álló eseménytér létrehozása
    end
    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
    y = lsim (system, szinusz);

    W = zeros(Pi,25000); % kezdeti mátrix nullázása

    for i=Pi:25000/2
        X=szinusz(i-Pi+1:i); % P darab minta kivétele
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*X'*e; % LMS algoritmus számolása
    end

    system = tf ([(1-(r-q)), 0, 0], [1, 0, 0, (r-q)], 1); % r q-val való csökkentése
    y = lsim (system, szinusz);

    for i=25000/2+1:25000
        X=szinusz(i-Pi+1:i); % P darab minta kivétele
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*X'*e; % LMS algoritmus számolása
    end

    fprintf('dimW=%dP eset együtthatói, eredeti',turn);
    W(:,25000/2)
    
    fprintf('dimW=%dP eset együtthatói, csökkentett r',turn);
    W(:,25000)
    
    figure(7+turn-1)
    [maxW, maxNum] = sort(abs(W(:,end-1)),'descend'); % eremények csökkenõ sorrendbe rendezése
    hold all
    for i=1:5
        plot(abs(W(maxNum(i),:)')); % 5 legnagyobb együttható kirajzolása
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

    lambda = 0.9; % feladat álltal megadott paraméterek
    v = 0.1;

    for i=2:M/2
        s = s+sin(i*2*pi*db+phase(i-1));
    end

    for i=1:99
        szinusz = [szinusz s]; % 25000 mintából álló eseménytér létrehozása
    end
    system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
    y = lsim (system, szinusz);

    W = zeros(Pi,25000);
    R = eye(Pi); % kezdeti mátrixok nullázása

    for i=Pi:25000/2
        X=szinusz(i-Pi+1:i); % P darab minta kivétele
        num = (R*X'*X*R);
        den = ((lambda/v)+X*R*X');
        R = (1/lambda)*(R-(num/den));
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*R*X'*e; % algoritmus számolása
    end

    system = tf ([(1-(r-q)), 0, 0], [1, 0, 0, (r-q)], 1); % r q-val való csökkentése
    y = lsim (system, szinusz);

    for i=25000/2+1:25000
        X=szinusz(i-Pi+1:i); % P darab minta kivétele
        num = (R*X'*X*R);
        den = ((lambda/v)+X*R*X');
        R = (1/lambda)*(R-(num/den));
        e=y(i)-X*W(:,i);
        W(:,i+1)=W(:,i)+2*LMS_mu*R*X'*e; % algoritmus számolása
    end
    
    fprintf('dimW=%dP eset együtthatói, eredeti',turn);
    W(:,25000/2)
    
    fprintf('dimW=%dP eset együtthatói, csökkentett r',turn);
    W(:,25000)
    
    figure(9+turn-1)
    [maxW, maxNum] = sort(abs(W(:,end-1)),'descend');
    hold all
    for i=1:5
        plot(abs(W(maxNum(i),:)')); % 5 legnagyobb együttható kirajzolása
    end
    legend(strcat('W',num2str(maxNum(1:5),-1)))
    hold off
end

%%%%%% 4. feladat %%%%%%%
P=6; % 6 együtthatót kell keresni

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
LMS_mu = 1/(1000*lambda); % bátorsági tényezõ meghatározása
szinusz = s;

for i=2:M/2
    s = s+sin(i*2*pi*db+phase(i-1));
end

for i=1:99
    szinusz = [szinusz s];
end
system = tf ([(1-r), 0, 0], [1, 0, 0, r], 1);
y = lsim (system, szinusz);

W = zeros(P,25000); % kezdeti mátrix nullázása

for i=P:25000-1
    X=[szinusz(i-1:i), y(i-3:i)']; % a megadott formula szerint 2 gerjesztõ minta,
                                   % valamint 4 kimeneti minta kiválasztása
    e=y(i+1)-X*W(:,i);
    W(:,i+1)=W(:,i)+2*LMS_mu*X'*e;
end

fprintf('Együtthatói:',turn);
W(:,end)

figure(11)
hold all
for i=1:6
	plot(W(i,:)');
end
legend('W0','W1','W2','W3','W4','W5',-1)
plot(W')
