%%
%
clc
clear
close all

x0       = 0000; 
y0       = 8000;
z0       = 0;
zsar     = 10000; % высота полета
R0       = sqrt(x0.^2 + y0^2 + zsar^2);

My    = 538;
gr       = 180 / pi;
Fprf = 700;
Tsyn = 0.77;
Tp   = 1/Fprf;
tn    = zeros(1, My);
R     = zeros(1, My);
ty   = (0 : My-1)*Tp;
Vsar  = 250;

y     = -10;
TetaQ = 90;
Rs     = zeros(1, My);
Rx     = zeros(1, My);

for ny = 1 : My
    tn(ny) = ty(ny) - Tsyn / 2;
    tnm = tn(ny) + y/Vsar;
    R(ny)  = sqrt(R0^2 + Vsar^2*tnm^2);
    TetaI  = atan(y/R0);
%     Rs(ny) =  R0 + Vsar^2*tn(ny)^2*sin(TetaQ/gr - TetaI)^2 / (2*R0);
    Rs(ny) = sqrt((x0-y)^2 + y0^2 + zsar^2 + 2*y*Vsar*tn(ny) + (Vsar*tn(ny))^2);
    Rx(ny) = sqrt((x0 - y - tn(ny)*Vsar)^2 + y0^2 + zsar^2);
    
end

figure
plot(1:My, Rs, 'x-r', 1:My, R, '.-b')
grid on