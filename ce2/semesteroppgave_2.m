clear;
defaultStream = RandStream.getDefaultStream();
reset(defaultStream);
%% Oppgave 1a
k = tan(pi.*[0.38,0.4,0.46,0.5]);
imphp = synN(1,0,k);
implp = synN(0,1,k)

figure(); stem(imphp); xlabel('n'); ylabel('h_{hp}(n)');
figure(); stem(implp); xlabel('n'); ylabel('h_{lp}(n)');

%% Oppgave 1b
% Kalkulerer og plotter korrelajsonene
[thp,Rhh]=xcorr(imphp,imphp);
[tlp,Rll]=xcorr(imphp,imphp);
[thl,Rhl]=xcorr(imphp,implp);

figure(); stem(Rhh,thp); xlabel('\tau (lag)'); ylabel('R_{hh}(\tau)');
figure(); stem(Rll,tlp); xlabel('\tau (lag)'); ylabel('R_{ll}(\tau)');
figure(); stem(Rhl,thl); xlabel('\tau (lag)'); ylabel('R_{hl}(\tau)');

%% Oppgave 1c
figure(); freqz(imphp);
figure(); freqz(implp);

%% Oppgave 1d
N = 10000

% Genererer gaussisk støy
En1g = randn(1,N);
En2g = randn(1,N);

% Modulerer
Yn1 = synN(En1g,En2g,k);
var(Yn1)

%% Oppgave 1e

% Genererer binær gaussisk støy
En1b = (randi(2,1,N)-1.5)*2;
En2b = (randi(2,1,N)-1.5)*2;

% Modulerer
Yn2 = synN(En1b,En2b,k);
var(Yn2)

%% Oppgave 1f

% Rekonstruerer signalene
[Xn1g,Xn2g] = anaN(Yn1,k);
[Xn1b,Xn2b] = anaN(Yn2,k);

% Finner maksimum differanse mellom innsignalene og de rekonstruerte
maxdiff_hpg = max( abs(En1g - Xn1g(5:N+4)) )
maxdiff_lpg = max( abs(En2g - Xn2g(5:N+4)) )
maxdiff_hpb = max( abs(En1b - Xn1b(5:N+4)) )
maxdiff_lpb = max( abs(En2b - Xn2b(5:N+4)) )

% [tx1g,Rx1g] = xcorr(En1g,Xn1g,25);
% [tx2g,Rx2g] = xcorr(En2g,Xn2g,25);
% [tx1b,Rx1b] = xcorr(En1b,Xn1b,25);
% [tx2b,Rx2b] = xcorr(En2b,Xn2b,25);
% 
% figure(); stem(Rx1g,tx1g); xlabel('\tau (lag)'); ylabel('R_{xx}(\tau)');
% figure(); stem(Rx2g,tx2b); xlabel('\tau (lag)'); ylabel('R_{xx}(\tau)');
% figure(); stem(Rx1b,tx1b); xlabel('\tau (lag)'); ylabel('R_{xx}(\tau)');
% figure(); stem(Rx2b,tx2b); xlabel('\tau (lag)'); ylabel('R_{xx}(\tau)');



%% Oppgave 2a

imp1 = synN(synN(1,0,k),synN(0,0,k),k);
imp2 = synN(synN(0,1,k),synN(0,0,k),k);
imp3 = synN(synN(0,0,k),synN(1,0,k),k);
imp4 = synN(synN(0,0,k),synN(0,1,k),k);

figure(); stem(imp1);
figure(); stem(imp2);
figure(); stem(imp3);
figure(); stem(imp4);

figure(); freqz(imp1);
figure(); freqz(imp2);
figure(); freqz(imp3);
figure(); freqz(imp4);

%% Oppgave 2b


