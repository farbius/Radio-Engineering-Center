%%
%
clc
clear
close all

c = 3e8;
dev  = 150e6;
dt   = 1/dev/8;
tau  = 6e-6;

R_0  = 14000;
tmin = 2*R_0/c - tau;
tmax = 2*R_0/c + tau ;
Mx   = 2*ceil((tmax - tmin)/2/dt);

R = 13900.00;
tx   = tmin + (0:Mx-1)*dt;

2*R/c/1e-6

td     = tx - 2*R/c; 
s_raw  = exp(1i*pi*dev/tau*(td.^2-td*tau)).*(td>=0 & td<=tau);


R3       = 14000.0;
td0      = tx - 2*R3/c;
h_range  = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
hF_range = fft(h_range);
fs_raw   = fft(s_raw);
fc_raw   = fs_raw.*conj(hF_range);
s_range  = fftshift(ifft(fc_raw));

figure
plot(tx.*c/2, abs(s_range), '.-b')
grid on