%% $$m(t)=2cos(200\pi t)$$
close all
clc
fmax = 100;
fs = 50*fmax;
ts = 1/fs;
T = 1/fmax;
t =-3*T:ts:3*T;
x = 2*cos (200*pi*t);
% Transformada
t=-50*T:ts:50*T;
M=0;
x = 2*cos (200*pi*t);
n=0;
w=-500*pi:500*pi;
for tt = t
   n =n+1;
   M= M +x(n)*exp(-j*w*tt);
end
ejex = [-.035 .035];
titulo = 'Señales en el dominio del tiempo y de frecuencia';
ejeX = '$t$';
ejeY = '$m(t)$';
simbo = '$m(t)= 2cos(200\pi t)$';
freq = w/(2*pi);
magn = abs (M);
ejex1 = [95 105];
ejex2 = [-250 250];
titulo1 = 'Delta de cosenos';
ejeX1 = '$\omega$';
ejeY1 = '$|M(\omega)|$';
simbo1= '$\frac{1}{2\pi}\int_{n=-\infty}^{\infty}F(\omega)e^{jn\omega_{0}t}d\omega$';
graficar2(t, x,ejex,titulo, ejeX, ejeY, simbo, freq, magn, ejex2,ejeX1, ejeY1, simbo1)
graficar(freq,magn,ejex1,titulo1, ejeX1, ejeY1, simbo1)
%% $$m(t) = 2cos(200pi*t) + 5sin(600\pi t)$$
clear 
close all
clc
fmax = 300;
fs = 50*fmax;
ts = 1/fs;
T = 1/fmax;
t =-3*T:ts:3*T;
x = 2*cos (200*pi*t) +5*sin(600*pi*t) ;
% Transformada
t=-50*T:ts:50*T;
M=0;
x = 2*cos (200*pi*t) +5*sin(600*pi*t);
n=0;
w=-1000*pi:1000*pi;
for tt = t
   n =n+1;
   M= M +x(n)*exp(-j*w*tt);
end
ejex = [-.035 .035];
titulo = 'Señales en el dominio del tiempo y de frecuencia';
ejeX = '$t$';
ejeY = '$m(t)$';
simbo = '$m(t)= 2cos(200\pi t) + 5sin(600\pi*t)$';
freq = w/(2*pi);
magn = abs (M);
ejex1 = [50 350];
a = 500;
ejex2 = [-a a];
titulo1 = 'Delta de cosenos';
ejeX1 = '$\omega$';
ejeY1 = '$|M(\omega)|$';
simbo1= '$\frac{1}{2\pi}\int_{n=-\infty}^{\infty}F(\omega)e^{jn\omega_{0}t}d\omega$';
graficar2(t, x,ejex,titulo, ejeX, ejeY, simbo, freq, magn, ejex2,ejeX1, ejeY1, simbo1)
graficar(freq,magn,ejex1,titulo1, ejeX1, ejeY1, simbo1)

%% $$m(t)=sin(\frac{20\pi*t}{20\pi*t})=sinc (20t)$$
clear 
close all
clc
fmax = 300;
fs = 50*fmax;
ts = 1/fs;
T = 1/fmax;
t =-50*T:ts:50*T;
x = sinc (20*t) ;
% Transformada
t=-50*T:ts:50*T;
M=0;
x = sinc (20*t) ;
n=0;
w=-100*pi:100*pi;
for tt = t
   n =n+1;
   M= M +x(n)*exp(-j*w*tt);
end
ejex = [-.2 .2];
titulo = 'Señales en el dominio del tiempo y de frecuencia';
ejeX = '$t$';
ejeY = '$m(t)$';
simbo = '$m(t)= sin(\frac{20\pi*t}{20\pi*t}) = sinc (20t)$';
freq = w/(2*pi);
magn = abs (M);
ejex1 = [-50 50];
a = 50;
ejex2 = [-a a];
titulo1 = 'Delta de cosenos';
ejeX1 = '$\omega$';
ejeY1 = '$|M(\omega)|$';
simbo1= '$\frac{1}{2\pi}\int_{n=-\infty}^{\infty}F(\omega)e^{jn\omega_{0}t}d\omega$';
graficar2(t, x,ejex,titulo, ejeX, ejeY, simbo, freq, magn, ejex2,ejeX1, ejeY1, simbo1)
graficar(freq,magn,ejex1,titulo1, ejeX1, ejeY1, simbo1)
%%
clear all
close all
clc
%%
function graficar(x,y,ejex,titulo, ejeX, ejeY, simbo)
    yy = .1*y;
    figure;
    plot (x, y,'k','LineWidth',1.5)
    title (titulo)
    xlabel (ejeX,'Interpreter','latex')
    ylabel (ejeY,'Interpreter','latex')
    axis([ejex(1) ejex(2) min(y-yy) max(y+yy)])
    legend(simbo,'Interpreter','latex','Location','Best')
    grid on
end
function graficar2(x, y,ejex,titulo, ejeX, ejeY, simbo, x1, y1, ejex1,ejeX1, ejeY1, simbo1)
    yy = .1*y;
    yy1 = .1*y1;
    figure;
    subplot(211)
    plot (x, y,'k','LineWidth',1.5)
    title (titulo)
    xlabel (ejeX,'Interpreter','latex')
    ylabel (ejeY,'Interpreter','latex')
    axis([ejex(1) ejex(2) min(y+yy) max(y+yy)])
    legend(simbo,'Interpreter','latex','Location','Best')
    grid on
    subplot(212)
    plot (x1, y1,'k','LineWidth',1.5)
    xlabel (ejeX1,'Interpreter','latex')
    ylabel (ejeY1,'Interpreter','latex')
    axis([ejex1(1) ejex1(2) min(y1-yy1) max(y1+yy1)])
    legend(simbo1,'Interpreter','latex','Location','Best')
    grid on
end