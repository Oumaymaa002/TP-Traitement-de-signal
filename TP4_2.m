clear all
close all
clc

[music,fs]=audioread("test.wav"); %lecture de l'audio

Ts= 1/fs; % Période d'échantillonage
N=length(music); % le nbr d'echantillons égal à la taille du vecteur music
t=0:Ts:(N-1)*Ts;
plot(t,music)
title('representation temporelle d un signal audio ')
xlabel('t')
xlim([1 1.2]);
grid on


fshift = (-N/2:N/2-1)*(fs/N); % le pas de discrétisation : fe/N
transfF=fft(music); % transformée de fourirer
plot(fshift,fftshift(abs(transfF)));
title('representation fréquentielle du signal audio ')
xlabel('f')
ylabel('Amplitude')

fc=4700; %fréquence de coupre, 4700 car l'atténuation est de 0.7
f=(0:N-1)*(fs/N);
k=1;

%Transmittance complexe
H=k./(1+j*(f/fc).^1000);

%Conception du filtre
H_filter = [H(1:floor(N/2)), flip(H(1:floor(N/2)))];

%Filtrage
f = transfF(1:end-1).*H_filter;

%Signal filtré
filtered_sign=ifft(f,"symmetric");
plot(fshift(1:end-1), fftshift(abs(fft(filtered_sign))));

