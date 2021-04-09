%% MODULACION AM 
% Con mu = .5
close all
clear 
clc
fs=50000;
ts=1/fs;
t0=0.15;
t=0:ts:t0;
m=1*(t>=0&t<t0/3)-2*(t>=t0/3&t<2*t0/3);

%% Señal $m(t)$
P=.5;
figure;
subplot (211);
plot (t,m,'k','LineWidth',1.5);
title ('Mensaje');
xlabel('$t$','Interpreter','latex');
ylabel('$m(t)$','Interpreter','latex');
axis([-.01 .16 -2.5 1.5])
grid on
legend('$m(t)$','Interpreter','latex','Location','Best')

subplot (212);
myTransformada (m,P);
title ('Transformada de M');
xlabel('$\omega$','Interpreter','latex');
ylabel('$M(\omega)$','Interpreter','latex');
axis([-500 500 -.005 .045])
legend('$|M(\omega)|$','Interpreter','latex','Location','Best')
grid on

%% Portadora $c(t)$
c=cos(1000*pi*t);
P=.8;
figure;
subplot (211);
plot (t,c,'k','LineWidth',1.1);
title ('Portadora');
xlabel('$t$','Interpreter','latex');
ylabel('$c(t)$','Interpreter','latex');
legend('$m(t)$','Interpreter','latex','Location','Best')
axis([-.005 .155 -1.5 1.5])
grid on

subplot (212);
myTransformada (c,P);
title ('Transformada de la portada');
xlabel('$\omega$','Interpreter','latex');
ylabel('$C(\omega)$','Interpreter','latex');
axis([-850 850 -.002 .022])
legend('$|M(\omega)|$','Interpreter','latex','Location','Best')
grid on
%% Modulación $Y_{AM}$
yam=m.*c;
mp=abs(min(yam));
mu=.5;
A=mp/mu;
yam=(A+m).*c;

P=1;
figure;
subplot (211);
plot (t,yam,'k','LineWidth',1.1);
title ('Modulacion $y_{am}(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$y_{am}(t)$','Interpreter','latex');
legend('$m(t)$','Interpreter','latex','Location','Best')
axis([-.005 .155 -6 6])
grid on

subplot (212);
myTransformada (yam,P);
title ('Transformada de Fourier $Y_{AM}(\omega)$','Interpreter','latex');
xlabel('$\omega$','Interpreter','latex');
ylabel('$Y_{AM}(\omega)$','Interpreter','latex');
legend('$|M(\omega)|$','Interpreter','latex','Location','Best')
axis([-1050 1050 -.01 .11])
grid on

%% Potencia de $Y_{AM}$ y $SNR = 10$ 
Pyam = potencia (yam);
SNR = 10;
Pr = Pyam*10^(-SNR/10);
sigma = sqrt (Pr);
ruido = sigma* randn (1, length(t));
Pruido = potencia (ruido);

figure;
subplot (211)
plot (t, ruido,'LineWidth',1.1)
title ('Ruido')
ylabel ('ruido')
xlabel ('t')
legend('Ruido','Location','Best')
axis([-.01 .16 min(ruido+.1*ruido) max(ruido+.1*ruido)])
grid on

subplot (212);
myTransformada (ruido,P);
title ('Transformada del ruido');
xlabel('$\omega$','Interpreter','latex');
ylabel('$R(\omega)$','Interpreter','latex');
legend('$R(\omega)$','Location','Best','Interpreter','latex')
grid on

mcanal = ruido + yam;

figure;
subplot (211)
plot (t, mcanal,'k','LineWidth',1.1)
title ('Señal con ruido')
ylabel ('$y_{am} + ruido$','Interpreter','latex')
xlabel ('$t$','Interpreter','latex')
axis([-.01 .16 min(mcanal+.1*mcanal) max(mcanal+.1*mcanal)])
legend('Señal con ruido','Location','Best')
grid on
subplot (212);
myTransformada (mcanal,P);
title ('Transformada de la señal con ruido');
xlabel('$\omega$','Interpreter','latex');
ylabel('$R(\omega)$','Interpreter','latex');
legend('$R(\omega)$','Location','Best','Interpreter','latex')
grid on
yam = mcanal;
%% Demodulacion coherente $r=y_{am}(t)*c(t)$
r=yam.*c;
P=2;
figure;
subplot (211);
plot (t,r,'k','LineWidth',1.1);
title ('Demodulacion coherente r');
xlabel('$t$','Interpreter','latex');
ylabel('$r(t)$','Interpreter','latex');
legend('$r(t)$','Location','Best','Interpreter','latex')
axis ([-.01 .16 min(r+.1*r) max(r+.1*r)])
grid on

subplot (212);
myTransformada (r,P);
title ('Transformada de r');
xlabel('$\omega$','Interpreter','latex');
ylabel('$R(\omega)$','Interpreter','latex');
legend('$R(\omega)$','Location','Best','Interpreter','latex')
grid on
%% Filtro pasabajas $h=800sinc(800t)$
h=800*sinc(800*t);
P=.8;
figure;
subplot (211);
plot (t,h,'k','LineWidth',1.5);
title ('Filtro h');
xlabel('$t$','Interpreter','latex');
ylabel('$h(t)$','Interpreter','latex');
legend('$h(t)$','Location','Best','Interpreter','latex')
axis ([-.01 .04 min(h+.1*h) max(h+.1*h)])
grid on

subplot (212);
myTransformada (h,P);
title ('Transformada de h');
xlabel('$\omega$','Interpreter','latex');
ylabel('$H(\omega)$','Interpreter','latex');
legend('$H(\omega)$','Location','Best','Interpreter','latex')
grid on
%% Aplicacion de filtro
m_rec=conv(r, h)*ts;
P=1.5;
tq=-length(t):length(t);
tq1(1,:)=tq(1:15001);
tq=tq1;
figure;
subplot(211)
plot (tq,m_rec,'k','LineWidth',1.1);
title ('Aplicacion del filtro');
xlabel('$t$','Interpreter','latex');
ylabel('$m_{rec}(t)$','Interpreter','latex');
axis ([-8000 2000 min(m_rec+.2*m_rec) max(m_rec+.1*m_rec)]);
legend('$m_{rec}(t))$','Location','Best','Interpreter','latex')
grid on

subplot (212);
M=myTransformada (m_rec,P);
title ('Transformada de h');
xlabel('$\omega$','Interpreter','latex');
ylabel('$M_{REC}(\omega)$','Interpreter','latex');
legend('$H(\omega)$','Location','Best','Interpreter','latex')
grid on

%% Detector de envolvente
figure;
r2=yam;
I=find(r2<0);
r2(I)=0;
plot(t,yam,'k--','LineWidth',1.1)
hold on
plot(t,r2,'r','LineWidth',1.2)
title ('Detector de envolvente');
xlabel('$t$','Interpreter','latex');
ylabel('$y_{am}(t)$','Interpreter','latex');
legend('$y_{am}(t)$','Envolvente','Location','Best','Interpreter','latex')
grid on
%% Modulacion $Y_{usb}(\omega)$ y $Y_{lsb}(\omega)$
close all
clear 
clc

fs=2500;
ts=1/fs;
t0=0.15;

t=0:ts:t0;
N=1000;
%% Señal $m(t)$
m=1*(t>=0&t<t0/3)-2*(t>=t0/3&t<2*t0/3);
figure
plot(t,m)

subplot (211)
plot (t,m,'k','LineWidth',1.1)
title ('Señal del mensaje');
xlabel('$t$','Interpreter','latex');
ylabel('$m(t)$','Interpreter','latex');
legend('$m(t)$','Location','Best','Interpreter','latex')
a = -.01;b = .16;
axis([a b min(m+.5*m) max(m+1)])
grid on
subplot(212)
MyTransformada2 (m,N);
title ('Transformada del Mensaje');
xlabel('$\omega$','Interpreter','latex');
ylabel('$M(\omega)$','Interpreter','latex');
legend('$M(\omega)$','Location','Best','Interpreter','latex')
a = 1500;
axis([-a a -.01 .12])
grid on

%% Portadora $c(t)$
c=cos(1000*pi*t);
N=10000;
figure;
subplot (211)
plot (t,c,'k','LineWidth',1.1)
title ('Señal portadora');
xlabel('$t$','Interpreter','latex');
ylabel('$c(t)$','Interpreter','latex');
legend('$c(t)$','Location','Best','Interpreter','latex')
a = -.01;b = .16;
axis([a b -1.1 1.1])
grid on
subplot(212)
MyTransformada2 (c,N);
title ('Transformada de portadora');
xlabel('$\omega$','Interpreter','latex');
ylabel('$C(\omega)$','Interpreter','latex');
legend('$C(\omega)$','Location','Best','Interpreter','latex')
a = 1300;
axis([-a a -.01 .085])
grid on

%% Mosulacion $y_{dsb-sc}(t)$
ydsb_sc = m.*c;
N=10000;
figure;
subplot (211)
plot (t,ydsb_sc,'k','LineWidth',1.1)
title ('Se\~nal modulada $y_{dsb-sc}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$y_{dsb-sc}(t)$','Interpreter','latex','FontSize',12);
grid on
legend('$y_{dsb-sc}(t)$','Location','Best','Interpreter','latex','FontSize',12)
a = -.01;b = .16;
axis([a b -2.5 2])
subplot(212)
MyTransformada2 (ydsb_sc,N);
title ('Transformada de $Y_{DSB-SC}$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex');
ylabel('$Y_{DSB-SC}(\omega)$','Interpreter','latex','FontSize',12);
a = 1300;
axis([-a a -.005 .065])
grid on
legend('$Y_{DSB-SC}(\omega)$','Location','Best','Interpreter','latex','FontSize',12)

%% $SNR = 10$ db
Pydsb_sc = potencia (ydsb_sc);
SNR = 10;
P=1;
Pr = Pydsb_sc*10^(-SNR/10);
sigma = sqrt (Pr);
ruido = sigma* randn (1, length(t));
Pruido = potencia (ruido);

figure;
subplot (211)
plot (t, ruido,'k','LineWidth',1.1)
title ('Ruido')
ylabel ('ruido')
xlabel ('$t$','Interpreter','latex')
legend('$r(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b -1 1])
subplot (212);
MyTransformada2 (ruido,N);
title ('Transformada del ruido');
xlabel('$\omega$','Interpreter','latex');
ylabel('$R(\omega)$','Interpreter','latex');
legend('$R(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.001 .01])
%% Canal = ruido + $y_{dsb-sc}(t)$
mcanal = ruido + ydsb_sc;

figure;
subplot (211)
plot (t, mcanal,'k','LineWidth',1.1)
title ('Señal con el ruido')
ylabel ('$f(t)$','Interpreter','latex')
xlabel ('$t$','Interpreter','latex')
legend('$f(t)+r(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b min(mcanal+.1*mcanal) max(mcanal+.1*mcanal)])
subplot (212);
MyTransformada2 (mcanal,N);
title ('Transformada de la señal con ruido');
xlabel('$\omega$','Interpreter','latex');
ylabel('$F(R(\omega))$','Interpreter','latex');
legend('$F(R(\omega))+R(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.001 .06])
%% Funcion de Hilbert
t=0:ts:t0;
c1=sin(1000*pi*t);
yh=imag(hilbert(m));

figure;
subplot (211)
plot (t,yh,'k','LineWidth',1.5)
title ('Funcion de hilbert $y_h$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$y_h(t)$','Interpreter','latex','FontSize',12);
legend('$y_h(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b min(yh+.1*yh) max(yh+.1*yh)])
subplot(212)
MyTransformada2 (yh,N);
title ('Transformada de $Y_H$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$Y_H(\omega)$','Interpreter','latex','FontSize',12);
legend('$Y_H(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.001 .12])
%% Filtro $h(t)$ pasa bajas
t=-3*t0:ts:3*t0;
h=600*sinc(300*t);
N=7000;
figure;
subplot (211)
plot (t,h,'k','LineWidth',1.5)
title ('Filtro h','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$h(t)$','Interpreter','latex','FontSize',12);
legend('$h(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.5;b = .5;
axis([a b min(h+.1*h) max(h+.1*h)])
subplot(212)
MyTransformada2 (h,N);
title ('Transformada del Filtro H','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$H(\omega)$','Interpreter','latex','FontSize',12);
legend('$H(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.1 2.5])

%% Señal $Y_{usb}(t)$
t=0:ts:t0;
yusb = m.*c - yh.*c1 + ruido;
N=7000;
figure;
subplot (211)
plot (t,yusb,'k','LineWidth',1.1)
title ('Se\~nal modulada $y_{usb}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$y_{usb}(t)$','Interpreter','latex','FontSize',12);
legend('$y_{usb}(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b min(yusb+.1*yusb) max(yusb+.1*yusb)])
subplot(212)
MyTransformada2 (yusb,N);
title ('Transformada de $Y_{USB}$','Interpreter','latex','FontSize',12);
xlabel('$\omega$');
ylabel('$Y_{USB}(\omega)$','Interpreter','latex','FontSize',12);
legend('$Y_{USB}(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.01 .12])
%% Recuperar Señal con $r(t)$ de $y_{usb}(t)$
t=0:ts:t0;
rt = yusb.*c;
N=7000;
figure;
subplot (211)
plot (t,rt,'k','LineWidth',1.1)
title ('Se\~nal modulada $r_t$ de $y_{usb}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$r_t(t)$','Interpreter','latex','FontSize',12);
legend('$r_t(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b min(rt+.1*rt) max(rt+.1*rt)])
subplot(212)
MyTransformada2 (rt,N);
title ('Transformada de $R_T$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$R_T(\omega)$','Interpreter','latex','FontSize',12);
legend('$R_T(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.01 .06])
%% Señal Recuperada de $y_{usb}(t)$
m_rec=conv(rt,h)*ts;
t1=-1313:1312;
N=1500;
figure;
subplot (211)
plot (t1,m_rec,'k','LineWidth',1.1)
title ('Se\~nal Recuperada $m_{rec}(t)$ de $y_{usb}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$m_{rec}(t)$','Interpreter','latex','FontSize',12);
legend('$m_{rec}(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -1400;b = 1400;
axis([a b min(m_rec+.1*m_rec) max(m_rec+.1*m_rec)])
subplot(212)
MyTransformada2 (m_rec,N);
title ('Transformada de la  Filtro Se\~nal Recuperada $M_{REC}(\omega)$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$M_{REC}(\omega)$','Interpreter','latex','FontSize',12);
legend('$M_{REC}(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.01 .12])
%% Señal  $y_{lsb}(t)$
t=0:ts:t0;
ylsb = m.*c + yh.*c1 + ruido;
N=7000;
figure;
subplot (211)
plot (t,ylsb,'k','LineWidth',1.1)
title ('Se\~nal modulada $y_{lsb}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$y_{lsb}(t)$','Interpreter','latex','FontSize',12);
legend('$y_{lsb}(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b min(ylsb+.1*ylsb) max(ylsb+.1*ylsb)])
subplot(212)
MyTransformada2 (ylsb,N);
title ('Transformada de $Y_{LSB}(t)$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$Y_{LSB}(\omega)$','Interpreter','latex','FontSize',12);
legend('$Y_{LSB}(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.01 .12])
%% Recuperar Señal con $r(t)$ de $y_{lsb}(t)$
rt = ylsb.*c;
N=7000;
figure;
subplot (211)
plot (t,rt,'k','LineWidth',1.1)
title ('Se\~nal recuperada $r_t$ de $y_{lsb}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$r_t(t)$','Interpreter','latex','FontSize',12);
legend('$r_t(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -.01;b = .16;
axis([a b min(rt+.1*rt) max(rt+.1*rt)])
subplot(212)
MyTransformada2 (rt,N);
title ('Transformada de $R_T(\omega)$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$R_T(\omega)$','Interpreter','latex','FontSize',12);
legend('$R_T(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.01 .06])
%% Señal Recuperada de $y_{lsb}(t)$
m_rec=conv(rt,h)*ts;
t1=-1313:1312;
N=1500;
figure;
subplot (211)
plot (t1,m_rec,'k','LineWidth',1.1)
title ('Se\~nal Recuperada $m_{rec}(t)$ de $y_{lsb}(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
ylabel('$m_{rec}(t)$','Interpreter','latex','FontSize',12);
legend('$m_{rec}(t)$','Location','Best','Interpreter','latex','FontSize',12)
grid on
a = -1400;b = 1400;
axis([a b min(m_rec+.1*m_rec) max(m_rec+.1*m_rec)])
subplot(212)
MyTransformada2 (m_rec,N);
title ('Transformada de la  Filtro Se\~nal Recuperada $M_{REC}$','Interpreter','latex','FontSize',12);
xlabel('$\omega$','Interpreter','latex','FontSize',12);
ylabel('$M_{REC}(\omega)$','Interpreter','latex','FontSize',12);
legend('$M_{REC}(\omega)$','Location','Best','Interpreter','latex','FontSize',12)
a = 1300;
axis([-a a -0.01 .12])
close all