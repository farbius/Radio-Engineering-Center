clc
clear
close all

dev = 50e6;
dt = 1/dev;
tau  = 4e-6;
tmin = 2e-6;
tmax = 10e-6;
Mx = 2*ceil((tmax - tmin)/2/dt);
tx = tmin + (0:Mx-1)*dt;
R = 600;
tz = 2*R/3e8;
td = tx - tz;
s = exp(1i*pi*dev/tau*(td.^2-td*tau)).*(td>=0 & td<=tau); %

td0 = tx - 6e-6;
s0 = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
fs0=fft(s0);

Corr = fft(s) .* conj(fs0);
smb = fftshift(ifft(Corr));

Corr0 = fft(s0) .* conj(fs0);
smb0 = fftshift(ifft(Corr0));

figure
% plot(tx.*1e6, real(s), 'x-r',tx.*1e6, real(s0), '.-b')
subplot(2,1,1)
plot(tx.*1e6, real(s), '.-b')
xlabel('\bf t_x,\rm ìêñ')
ylabel('\rm real(\bf U_0_i_k)')
ylim([-1.2 1.2])
xlim([tmin.*1e6 tmax.*1e6])
% hold on
% plot(tx.*1e6, imag(s), 'x-g',tx.*1e6, imag(s0), '.-y')
grid on
subplot(2,1,2)
plot(tx.*1e6, imag(s), '.-b')
xlabel('\bf t_x, \rm ìêñ')
ylabel('\rm imag(\bf U_0_i_k)')
ylim([-1.2 1.2])
xlim([tmin.*1e6 tmax.*1e6])
grid on

%% Ğàñ÷åò ÷àñòîòû Äîïëåğà
gr   = 180 / pi;
Vsar = 250; 
zsar = 10000;
Tsyn = 1.0;
Tp   = 2e-3;
x0   = 4000;
y0   = 20000;
z0   = 0;
Lam  = 0.03;

Teta = 90 - atan(x0/y0).*gr

My   = 2*round(.5*(Tsyn/Tp));
xsar = (-My/2 : My/2)*Tp*Vsar;
ty   = (-My/2 : My/2)*Tp;

R    = sqrt((x0 - xsar).^2 + y0^2 + zsar^2);

R0  = sqrt(x0.^2 + y0^2 + zsar^2);
La   = 2; 
Fdna = sinc(La*atan(Vsar.*ty./R0)/Lam);


SDop = exp(j*4*pi*R./Lam);
sF   = abs(fftshift(fft(SDop)));
F    = (-My/2 : My/2)/Tsyn;

figure
subplot(2,1,1)
plot(ty, real(SDop), '.-b', ty, imag(SDop), '.-r', ty, Fdna, '.-k')
ylim([-1.2 1.2])
xlabel('\bf T_{syn}, ìñ')
ylabel('\bf s_\Omega_i{_k}')
legend('\rm real(\bf s_\Omega_i{_k})', '\rm imag(\bf s_\Omega_i{_k})', '\bf Fdna')
grid on
subplot(2,1,2)
plot(F, sF, '.-b')
xlabel('\bf f, \bf Ãö')
ylabel('\bf |sF_\Omega_i_k|')
grid on