clc
clear
close all

%% –асчет частоты ƒоплера
gr   = 180 / pi;
Vsar = 250; 
zsar = 10000;
Tsyn = 1.0;
Tp   = 2e-3;
x0   = 4000;
y0   = 20000;
z0   = 0;
Lam  = 0.03;
%% блок координат наземных целей
x1   = -0;  % м, смещение X относительно центра участка картографировани€
y1   = -50; % м, смещение Y относительно центра
%%

Teta = 90 - atan(x0/y0).*gr;

My   = 2*round(.5*(Tsyn/Tp));
xsar = (-My/2 : My/2)*Tp*Vsar;
ty   = (-My/2 : My/2)*Tp;

R    = [sqrt((x0 - xsar).^2 + y0^2 + zsar^2);
        sqrt((x0 - xsar + x1).^2 + (y0 + y1)^2 + zsar^2)];

R0  = sqrt(x0.^2 + y0^2 + zsar^2);
La   = 2; 
Fdna = [sinc(La*atan(Vsar.*ty./R0)/Lam);
        sinc(La*atan(Vsar.*(ty-y1/Vsar)./R0)/Lam)];

SDop = zeros(size(R));
sF   = zeros(size(R));

for i = 1:size(R, 1)
SDop(i, :) = Fdna(i, :).*exp(1i*4*pi*R(i, :)./Lam);
  sF(i, :) = abs(fftshift(fft(SDop(i, :))));
end
F    = (-My/2 : My/2)/Tsyn;

figure
subplot(2,1,1)
plot(ty, real(SDop(1, :)), '.-b', ty, real(SDop(2,:)), '.-r')
ylim([-1.2 1.2])
xlabel('\bf T_{syn}, мс')
ylabel('\bf real(s_\Omega_i{_k})')
legend('\rm real(\bf s_\Omega_1{_k})', '\rm real(\bf s_\Omega_2{_k})')
grid on
subplot(2,1,2)
plot(F, sF(1, :), '.-b', F, sF(2, :), '.-r')
xlabel('\bf f, \bf √ц')
ylabel('\bf |sF_\Omega_i_k|')
legend('\rm \bf |sF_\Omega_1{_k}|', '\rm \bf |sF_\Omega_2{_k}|')
grid on