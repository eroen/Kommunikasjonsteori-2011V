%% Computer Assignment 2
% TTT4115 Kommunikasjonsteori
%
% By Erik Moen, Geir Huth and Qingrui Zhou.

%% Oppgave 1a
% Impulse responses for the high- and lowpass channels.

clear;
close all;
defaultStream = RandStream.getDefaultStream();
reset(defaultStream);

k = tan(pi.*[0.38,0.4,0.46,0.5]);
imphp = synN(1,0,k);
implp = synN(0,1,k);

figure(); stem(imphp); xlabel('n'); ylabel('h_{hp}(n)'); 
  title('Impulse response for the highpass channel');
figure(); stem(implp); xlabel('n'); ylabel('h_{lp}(n)'); 
  title('Impulse response for the lowpass channel');

%% Oppgave 1b
% We show the impulses to be orthogonal at the required offsets, first visually.
[Rhh,thp]=xcorr(imphp,imphp);
[Rll,tlp]=xcorr(imphp,imphp);
[Rhl,thl]=xcorr(imphp,implp);

figure(); stem(thp,Rhh); xlabel('\tau (lag)'); ylabel('R_{hh}(\tau)'); 
title('Autocorrelation for the highpass channel impulse response');
figure(); stem(tlp,Rll); xlabel('\tau (lag)'); ylabel('R_{ll}(\tau)'); 
title('Autocorrelation for the lowpass channel impulse response');
figure(); stem(thl,Rhl); xlabel('\tau (lag)'); ylabel('R_{hl}(\tau)'); 
title('Cross correlation between the high- and lowpass channel impulse responses');
snapnow;


%%
% The maximal unwanted correlations are shown below. As we can see, there are only residual roundoff errors.
center = find(thp == 0);
idx=sort([center+(2:2:(center-1)), center-(2:2:(center-1))]);
max(abs(Rhh(idx)))
center = find(tlp == 0);
idx=sort([center+(2:2:(center-1)), center-(2:2:(center-1))]);
max(abs(Rll(idx)))
center = find(thl == 0);
idx=sort([center+(2:2:(center-1)), center-(2:2:(center-1)), center]);
max(abs(Rhl(idx)))


%% Oppgave 1c
% The frequency responses for the two channels are shown.
figure(); freqz(imphp);
title('Frequency response for highpass channel');
figure(); freqz(implp);
title('Frequency response for lowpass channel');

%% Oppgave 1d
% The variance over the channel is computed using two gaussian sources.
N = 10000;

% Genererer gaussisk støy
En1g = randn(1,N);
En2g = randn(1,N);

% Modulerer
Yn1 = synN(En1g,En2g,k);
var(Yn1)

%%
% We observe the variance to be the same as the source signals, which makes sense as the rate is doubled.

%% Oppgave 1e
% The variance over the channel is computed using two binary sources.

% Genererer binær gaussisk støy
En1b = (randi(2,1,N)-1.5)*2;
En2b = (randi(2,1,N)-1.5)*2;

% Modulerer
Yn2 = synN(En1b,En2b,k);
var(Yn2)

%%
% We observe the variance to be the same as the source signals, which makes sense as the rate is doubled.

%% Oppgave 1f
% We will show that the input signals can be reconstructed perfectly by subtracting the reconstructed output from the input.
%
% Note that the modulation and demodulation creates an offset of 4 samples.

% Rekonstruerer signalene
[Xn1g,Xn2g] = anaN(Yn1,k);
[Xn1b,Xn2b] = anaN(Yn2,k);

% Finner maksimum differanse mellom innsignalene og de rekonstruerte
maxdiff_hpg = max( abs(En1g - Xn1g(5:N+4)) )
maxdiff_lpg = max( abs(En2g - Xn2g(5:N+4)) )
maxdiff_hpb = max( abs(En1b - Xn1b(5:N+4)) )
maxdiff_lpb = max( abs(En2b - Xn2b(5:N+4)) )

%%
% We attribute the differences to rounding errors.

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
% The impulse and frequency responses for the channels in a 4-channel transmitter are drawn.
close all;

imp1 = synN( synN(1,0,k), synN(0,0,k) ,k)';
imp2 = synN( synN(0,1,k), synN(0,0,k) ,k)';
imp3 = synN( synN(0,0,k), synN(1,0,k) ,k)';
imp4 = synN( synN(0,0,k), synN(0,1,k) ,k)';
imp = [imp1, imp2, imp3, imp4];

figure(); stem(imp1); xlabel('n'); ylabel('h_1(n)');
title('Impulse response for channel 1');
figure(); stem(imp2); xlabel('n'); ylabel('h_2(n)');
title('Impulse response for channel 2');
figure(); stem(imp3); xlabel('n'); ylabel('h_3(n)');
title('Impulse response for channel 3');
figure(); stem(imp4); xlabel('n'); ylabel('h_4(n)');
title('Impulse response for channel 4');

figure(); freqz(imp1);
title('Frequency response for channel 1');
figure(); freqz(imp2); 
title('Frequency response for channel 2');
figure(); freqz(imp3);
title('Frequency response for channel 3');
figure(); freqz(imp4);
title('Frequency response for channel 4');

%% Oppgave 2b
% The orthogonality is tested.
%
% First, we  test the orthogonality of each channel against itself. As earlier, the requirement is that the channel is orthogonal for even offsets greater than 0.

[R, tau] = xcorr(imp);
center=find(tau == 0);
autocol=5*(0:3)+1;

autotestrow=sort([center+(4:4:center-1), center-(4:4:center-1)]);

max(max(abs(R(autotestrow, autocol))))
%%
% We can see the maximal correlation is very small, an we attribute it to numerical error.
%
% Next, we look at orthogonality between channels. We demand that the channels are orthogonal for offsets that are multiples of 4.

crosscol=setdiff(1:16, autocol);
crosstestrow=sort([center+(4:4:center-1), center-(4:4:center-1), center]);

max(max(abs(R(crosstestrow, crosscol))))
%%
% We observe this number to also be very small, and attribute the value to rounding error.

%% Problem 2 c)
% We shall verify the integrity of a decoded signal.

defaultStream = RandStream.getDefaultStream();
reset(defaultStream);
L=1e4;
x=randn(4, L);
y=synN( synN(x(1,:), x(2,:), k), synN(x(3,:), x(4,:), k), k);
%[zprime]=zeros(2, 20016);
[zprime(1,:), zprime(2,:)]=anaN(y, k);
[z(1,:), z(2,:)]=anaN(zprime(1,:), k);
[z(3,:), z(4,:)]=anaN(zprime(2,:), k);

%%
% First, we must find the offset.
offset=find(xcorr(x(1,:), z(1,:))>=L/2)-L

%%
% To test for corruption, we will subtract the offset output from the input.
max(max(x - z(:, 1+offset:L+offset)))

%%
% We observe only a very small number, and consider the signal uncorrupted.

%% Problem 2 d)
% In order to introduce noise to match a SNR given in dB, we must first convert the value to a number using the |db2pow| function. We can then calculate the variance (power) of the transmitted signal to obtain the correct noise variance or power.

x=(randi(2, 4, L) * 2) - 3;
y=synN( synN(x(1,:), x(2,:), k), synN(x(3,:), x(4,:), k), k);
n0=randn(size(y));
S=var(y);

SNR=db2pow(4);
% We might have to divide this by 4.
N=S/SNR;
n=(N)*n0;
ynoisy=y+n;

%% Problem 2 e)
% We will determine the performance of the system by testing it for bit error rate with different SNRs.
snroffset=-30;
snrscale=0.5;
ival= 1:80;
for i = ival
	i
    SNR=db2pow(i*snrscale+snroffset);
    N=S/SNR;
    n=N*n0;
    yn=y+n;
    %[zprime]=zeros(2, 20016);
    [zprime(1,:), zprime(2,:)]=anaN(yn, k);
    [z(1,:), z(2,:)]=anaN(zprime(1,:), k);
    [z(3,:), z(4,:)]=anaN(zprime(2,:), k);

    zregen=zeros(size(z))-1;
    idx=find(z>0);
    %Should address the issue of exact 0.
    zregen(idx)=1;

    ierr=find(not(x == zregen(:, offset+1:L+offset)));
    nerr=length(ierr);
    BER(i)=nerr/(L*4);
    BER(i);
end

semilogy(ival*snrscale+snroffset, BER);



snapnow;
close all;
