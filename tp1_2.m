clear all
close all
clc

s='bluewhale.wav';
[son,Fe]=audioread(s);
sound(son,Fe);

son1 = son(2.45e4: 3.10e4);
sound(son,Fe);

plot(son1);
xlabel('t');
ylabel('son');
title('Représentation temporelle du son roqual bleu')
grid on

N = size(son1);

%question 3:

    fourier = fft(son1); % transformation de fourier rapide sur le chant 
    
    Densite_spectrale = abs(fourier).^2/N; %Densité spectrale du Chant
    
    f = (0:floor(N/2))*(Fe/N)/10;
    plot(f,Densite_spectrale(1:floor(N/2)+1));
    legend("Densité spectrale du chant");
    xlabel("Fréquence (Hz)");
    ylabel("Densité spectrale en puissance");