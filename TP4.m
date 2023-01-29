clear all
close all
clc

Te=0.0005;
fe = 1/Te;
t = 0:Te:5;

f1=500;
f2=400;
f3=50;
pi=3.14159265359;

%Représentation temporelle du signal
N=length(t);
x=sin(2*pi*f1*t)+sin(2*pi*f2*t)+sin(2*pi*f3*t);
plot(t,x);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal');

%Représentation fréquentielle du signal
fshift=(-N/2:N/2-1)*fe/N;
y=fft(x);
plot(fshift,fftshift(abs(y)));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal');



f = (0:N-1)*(fe/N);
K=1;
w=2*pi*f;
wc=50;



% semilogx(f,abs(H));
% xlabel('f');
% ylabel('Amplitude');
% title('Représentation du le module de la transmittance complexe');

wc=50;
wc1=500;
wc2=1000;

H=(K*j*w/wc)./(1+j*w/wc);
H1=(K*j*w/wc1)./(1+j*w/wc1);
H2=(K*j*w/wc2)./(1+j*w/wc2);

G=20*log(abs(H));
G1=20*log(abs(H1));
G2=20*log(abs(H2));

 subplot(211)
 semilogx(f,abs(H));

 subplot(212)
 semilogx(f,G,f,G1,f,G2) % Représentation du gain

 

%Filtrage du signal

filtre1 = H.*y;
filtre2 = H1.*y;
filtre3 = H2.*y;

x1=ifft(filtre1,"symmetric");
x2=ifft(filtre2,"symmetric");
x3=ifft(filtre3,"symmetric");

plot(fshift, fftshift(abs(fft(x3))));

subplot(411)
plot(t,x)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal initial');

subplot(412)
plot(t,x1)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal pour wc = 50');

subplot(413)
plot(t,x2)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal wc = 500');

subplot(414)
plot(t,x3)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal wc = 1000');


