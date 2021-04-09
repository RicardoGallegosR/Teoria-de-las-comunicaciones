%% Señal cuadrada
clear all
close all
clc
fc=500;
fs=50*fc;
ts=1/fs;
t0=0.15;
N=100000;
t=0:ts:t0;
f=linspace(-fs/2,fs/2, N);

%% Señal $m(t)$
m=1*(t>=0 & t<t0/3) - 2*(t>=t0/3 & t<2*t0/3);

figure
subplot (211)
plot(t,m,'k','LineWidth',1.5)
title ('Se\~nal $m(t)$','Interpreter','latex')
ylabel ('m(t)','Interpreter','latex')
xlabel ('t','Interpreter','latex')
axis([0, t0 -2.5 1.5])
grid on 
legend('$m(t)$','Location','Best','Interpreter','latex','FontSize',12)
M=fftshift(fft(m,N))*ts;
subplot (212)
k=abs(M);
M=k;
plot(f, M,'k','LineWidth',1.5)
axis([-200 200 1.1*min(M) 1.1*max(M)])
title ('Transformada de $|M(\omega)|$','Interpreter','latex')
ylabel ('$|M(\omega)|$','Interpreter','latex')
xlabel ('$\omega$','Interpreter','latex')
grid on
legend('$|M(\omega)|$','Location','Best','Interpreter','latex','FontSize',12)
%% Modulación en FM
Im=cumsum(m)*ts;
kf=100*pi;
A=1;
yfm=A*cos(2*pi*fc*t + kf*Im);
figure
subplot (211)
plot(t,yfm,'k','LineWidth',1.5)
title ('Modulacion en FM','Interpreter','latex')
ylabel ('m(t)','Interpreter','latex')
xlabel ('t','Interpreter','latex')
axis ([-.01 .16 -1.5 1.5])
grid on 
legend('$FM(t)$','Location','Best','Interpreter','latex','FontSize',12)

% derivada
dyfm=yfm(1);
n=1;
for tt=0:ts:t0-ts
    n=n+1;
    dyfm(n)=(yfm(n)-yfm(n-1))/ts;
end
th=-3*t0:ts:3*t0;
fcor=200;
h=pi*fcor*sinc(fcor*th);
subplot (212)
plot(t,dyfm,'k','LineWidth',1.5)
title ('Diferencial de FM','Interpreter','latex')
H=fftshift(fft(h,N))*ts;
grid on 
legend('$dyFM$','Location','Best','Interpreter','latex','FontSize',12)
%% Filtro pasabanda
figure
subplot (211)
plot(th, h,'k','LineWidth',1.5)
axis([-.05 .05 1.1*min(h) 1.1*max(h)])
title ('Filtro','Interpreter','latex')
ylabel ('Sa(t)','Interpreter','latex')
xlabel ('t','Interpreter','latex')
legend('$Sa(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on 
subplot (212)
k=abs (H);
H=k;
plot(f, abs(H),'k','LineWidth',1.5)
axis([-200 200 1.1*min(H) 1.1*max(H)])
title ('Transformada del Filtro $|H(\omega)|$','Interpreter','latex')
ylabel ('$|H(\omega)|$','Interpreter','latex')
xlabel ('$\omega$','Interpreter','latex')
grid on 
legend('$|H(\omega)|$','Location','Best','Interpreter','latex','FontSize',12)

%% Detector de envolvente
r=dyfm;
I=find(r<0);
r(I)=0;

figure
subplot (211)
plot(t,r,'k','LineWidth',1.5)
title ('Detector de envolvente','Interpreter','latex')
ylabel ('$|r(t)|$','Interpreter','latex')
xlabel ('t','Interpreter','latex')
axis([-.01 .16 -500 4000])
grid on 
legend('$|r(t)|$','Location','Best','Interpreter','latex','FontSize',12)

rr=conv(r,h)*ts;

subplot (212)
trr=th(1)+t(1):ts:th(end)+t(end);
plot(trr,rr,'k','LineWidth',1.5)
title ('Se\~nal recuperada','Interpreter','latex')
axis([-.1 .25 1.5*min(rr) 1.1*max(rr)])
ylabel ('$m_{rec}(t)$','Interpreter','latex')
xlabel ('t','Interpreter','latex')
legend('$m_{rec}(t)|$','Location','Best','Interpreter','latex','FontSize',12)

grid on 

%% Comparacion de señales original VS recuperada
mrec=(rr-A*2*pi*fc)/(A*kf);
figure
plot(trr,mrec,'r','LineWidth',1.5)
axis([0, t0 -2.4 1.1*max(mrec)])
hold on
plot(t,m,'k','LineWidth',1.5)
grid on
title ('Comparacion entre se\~nales','Interpreter','latex','FontSize',12)
ylabel ('$f(t)$','Interpreter','latex','FontSize',12)
xlabel ('t','Interpreter','latex','FontSize',12)
legend('$m_{rec}(t)|$','m(t)','Location','Best','Interpreter','latex','FontSize',12)
%%
close all