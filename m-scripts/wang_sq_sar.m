%%
%
clc
clear
close all

fc   = 10e10; % Hz
c    = 3e8;
Lam  = c/fc;

Tr   = 6.033e-6; % LFM time duration
Kr   = 4e12;     % Hz/s
dev  = Tr*Kr;    % LFM

R_0  = 7500;     % m
fs   = 30e6;     % sample rate (range)
Vsar = 200;      % m/s

Ka   = -2*Vsar/(R_0*Lam);
fprf = 500;      % Fprf = 500 Hz
La   = 1;        % antenna length

%% computed data
Rs   = c/fs;     % sample spacing
As   = Vsar/fprf;% azimuth spacing
TetaH= Lam/La;

Ls   = R_0*TetaH;% length aperture
Fdop = 2*Vsar/La;

Na   = ceil(Ls/As);
Nr   = ceil(fs/Tr);



