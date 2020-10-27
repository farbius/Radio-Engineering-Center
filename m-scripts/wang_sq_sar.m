%%
%
clc
clear
% close all

fc   = 10e9; % Hz
c    = 3e8;
gr   = 180 / pi;

Lam  = c/fc;
fprintf(">> длина волны %2.2f м \n", Lam);

Tr   = 6.033e-6; % LFM time duration
Kr   = 4e12;     % Hz/s
dev  = Tr*Kr;    % LFM
fprintf(">> девиация ЛЧМ %2.2f МГц \n", dev./1e6);

R_0  = 7500;     % m
fs   = 30e6;     % sample rate (range)
Vsar = 300;      % m/s

Ka   = -2*Vsar/(R_0*Lam);
fprf = 2000;     %
La   = 1;        % antenna length

TetaQ= 10;
fprintf(">> угол наблюдения %2.2f град \n", TetaQ);
%% computed data
Rs   = c/fs;     % sample spacing
As   = Vsar/fprf;% azimuth spacing
TetaH= Lam/La*gr;
fprintf(">> ширина гл луча ДНА %2.2f град \n", TetaH);


Ls   = R_0*TetaH/gr;% length aperture
Fdop = 2*Vsar/La;

Na   = ceil(Ls/As);
Nr   = ceil(fs/Tr);

fDL  = 2*Vsar/Lam*sin((TetaQ-.5*TetaH)/gr);
fprintf(">> нижний Доплер %8.2f Гц \n", fDL);
fDU  = 2*Vsar/Lam*sin((TetaQ+.5*TetaH)/gr);
fprintf(">> верхний Доплер %8.2f Гц \n", fDU);
fD   = fDU - fDL;
fDC  = fDL + fD/2;
fprintf(">> средний Доплер %8.2f Гц \n", fDC);
fprintf(">> ширина спектра тр сигнала %8.2f Гц \n", fD);
%% положение РСА
xU   = -R_0*tan((TetaQ+.5*TetaH)/gr);
xC   = -R_0*tan((TetaQ)/gr);
xL   = -R_0*tan((TetaQ-.5*TetaH)/gr);
Xsar = linspace(xU, xL, Na);
R    = sqrt(Xsar.^2 + R_0^2);
RC   = sqrt(xC.^2   + R_0^2);
RL   = sqrt(xL.^2   + R_0^2);
RU   = sqrt(xU.^2   + R_0^2);

Teta = linspace(TetaQ+.5*TetaH, TetaQ-.5*TetaH, Na);
f_dop= 2*Vsar/Lam*sin(Teta./gr);

figure(1)
plot(R, f_dop, '.-b', RL, fDL, 'or', RU, fDU, 'or', RC, fDC, 'or')
txt = ['Teta_H: ' num2str(TetaQ) ' град'];
text(RC+1.0,fDC,txt)
title('Изменение Доплера на интервале синтезирования')
xlabel('Дальность, м')
ylabel('Доплер, Гц')
hold on
grid on



