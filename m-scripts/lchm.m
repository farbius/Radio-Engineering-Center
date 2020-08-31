%%
%
clc
clear
close all

Fd = 100e6;
T0 = 10e-6;
dF = 10e6;
w0 = 2*pi*0e6;

N = round(T0*Fd);
% t = (-N/2:N/2)/Fd;
t = (0:N-1)/Fd;

Kr = pi*dF / T0;


s = cos(w0.*t + Kr.*t.^2/2);

figure
plot(t./1e-6,s, '.-')
xlabel('t, usec')
grid on