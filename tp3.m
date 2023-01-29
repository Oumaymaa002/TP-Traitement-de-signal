clear all
close all
clc

% Représentation temporelle

load('ecg.mat');
Fe=500;
Te=1/Fe;
N=length(ecg);
t=0:Te:(N-1)*Te;
plot(t,ecg)
grid on
title(" représentation temporelle de l’activation électrique du cœur")
xlabel("t")
ylabel("ECG")
xlim([0.5 1.5]);

% Représentation Fréquentielle

f=(0:N-1)*(Fe/N);
fshift=(-N/2:N/2-1)*Fe/N;
y = fft(ecg);
plot(fshift,fftshift(abs(y)));
grid on
xlabel('f');
ylabel('Amplitude')
title('Représentation temporelle de l’activation électrique du cœur');

% Conception du fitltre pass_haut

pass_haut = ones(size(ecg));
fc1=0.5;
index_fc = ceil((fc1*N)/Fe);
pass_haut(1:index_fc)=0;
pass_haut(N-index_fc+1:N)=0;
plot(f,pass_haut,"linewidth",1.5)
xlabel('f');
ylabel('Amplitude')
title('Conception du filtre');


% Filtrage
ecg1_freq = pass_haut.*y;
ecg1=ifft(ecg1_freq,"symmetric");
plot(t,ecg);
hold on
plot(t,ecg1+3);
hold on
plot(t,ecg-ecg1+1.5);
grid on
xlabel('t');
ylabel('signal')
title('Représentation du signal avant et après le filtre passe-haut');


% Conception du fitltre pass_notch

pass_notch=ones(size(ecg));
fc2=50;
index_fc2= ceil((fc2*N)/Fe)+1;
pass_notch(index_fc2)= 0;
pass_notch(N-index_fc2+1) = 0;

%Filtrage 2 

ecg2_freq = pass_notch.*fft(ecg1); 
ecg2 = ifft(ecg2_freq,"symmetric");

subplot(211)
%plot(t,ecg2)
plot(t,ecg-ecg1)
grid on
xlabel('t');
ylabel('signal')
title('Partie supprimée après application du filtre passe-haut');
subplot(212)
plot(t,ecg-ecg2)
grid on
xlabel('t');
ylabel('signal')
title('Partie supprimée après application du filtre notch');


plot(t,ecg);
hold on
plot(t,ecg2);
hold on
plot(t,ecg-ecg2);
grid on
xlabel('t');
ylabel('signal')
title('Partie supprimée après application du filtre notch');

% Conception du filtre pass_bas 

pass_bas=zeros(size(ecg2));
fc3 = 15;
%fc3 =10;
%fc3 = 20;
index_fc3 = ceil((fc3*N)/Fe);
pass_bas(1:index_fc3)= 1;
pass_bas(N-index_fc3+1:N) = 1;
plot(f,pass_bas,"linewidth",1.5)
xlabel('f');
ylabel('Amplitude')
title('Filtre pass-bas');

% Filtrage 3 

ecg3_freq = pass_bas.*fft(ecg2); 
ecg3 = ifft(ecg3_freq,"symmetric");

subplot(311)
plot(t,ecg)
xlabel('t');
ylabel('signal')
title('Signal ecg original');


subplot(312)
plot(t,ecg3)
xlabel('t');
ylabel('signal')
title('Signal ecg après application du filtre passe-bas');

subplot(313)
plot(t,ecg-ecg3)
xlabel('t');
ylabel('signal')
title('Partie supprimée après application du filtre passe-bas');

%Définir l'intervalle de recherche pour la fréquence cardiaque
min_hr = 40; % battements par minute
max_hr = 220; % battements par minute

% Calculer l'autocorrélation du signal ECG
[acf,lags] = xcorr(ecg3,ecg3);

% Trouver la fréquence cardiaque en se basant sur l'autocorrélation
[max_corr, max_index] = max(acf);
heart_rate = 60*Fe/(lags(max_index))

% Vérifier si la fréquence cardiaque est dans l'intervalle de recherche
if heart_rate > min_hr && heart_rate < max_hr
    disp(['Fréquence cardiaque : ', num2str(heart_rate), ' battements par minute']);
else
    disp('Fréquence cardiaque non détectée');
end
