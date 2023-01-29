clear all
close all
clc 

fe=10000;
Te=1/fe;
N=5000;
t = 0:Te:(N-1)*Te;
pi=3.14159265359;
%Représentation temporelle du signal

x = 1.2*cos(2*3.141592*440*t+1.2)+3*cos(2*3.141592*550*t)+0.6*cos(2*3.141592*2500*t);
x2 = 1.2*cos(2*pi*440*t+1.2)+3*cos(2*pi*550*t);
%%subplot(2,1,1)
%plot(t,x)
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal');
sound(x)

% Spectre du signal

f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*(fe/N);
y=fft(x);

%plot(fshift,fftshift(abs(y)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal');

% Ajout du bruit

xnoise = x+2*randn(size(x));
plot(t,xnoise)
plot(t,xnoise)
xlabel('t');
ylabel('xnoise');
title('Représentation temporelle du signal bruité');
%sound(xnoise)

%Spectre du signal bruité
f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*(fe/N);
ynoise=fft(xnoise);
plot(fshift,fftshift(abs(ynoise)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal bruité');

xnoise2 = x+10*randn(size(x));
plot(t,xnoise2)
xlabel('t');
ylabel('xnoise2');
title('Représentation temporelle du signal bruité');
%sound(xnoise)

%Spectre du signal bruité
f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*(fe/N);
ynoise2=fft(xnoise2);
plot(fshift,fftshift(abs(ynoise2)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal bruité');

subplot(321)
plot(t,x)
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal');

subplot(322)
plot(fshift,fftshift(abs(y)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal');

subplot(323)
plot(t,xnoise)
xlabel('t');
ylabel('xnoise');
title('Représentation temporelle du signal bruité');


subplot(324)
plot(fshift,fftshift(abs(ynoise)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal bruité');

subplot(325)
plot(t,xnoise2)
xlabel('t');
ylabel('xnoise2');
title('Représentation temporelle du signal bruité');

subplot(326)
plot(fshift,fftshift(abs(ynoise2)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal bruité');

%Conception du filtre

fc=2500
pass_bas=zeros(size(x));
index_fc = ceil((fc*N)/fe);
pass_bas(1:index_fc)= 1;
pass_bas(N-index_fc+1:N) = 1;
xlabel('f');
ylabel('Amplitude')
title('Filtre pass-bas');
plot(f,pass_bas);

% Filtrage 

sign_freq = pass_bas.*y; 
%plot(f,sign_freq);
filtered_sign = ifft(sign_freq,"symmetric");
plot(fshift,fftshift(abs(fft(filtered_sign))));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal filtré');


subplot(211)
plot(t,x)
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal iniltial');
subplot(212)
plot(t,filtered_sign)
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal filtré');
sound(filtered_sign)