%%
%

clc
clear
close all

fprintf("Start SAR Model \n");
%% SAR parameters
gr   = 180 / pi;
c    = 3e8;
Vsar = 250; 
zsar = 10000;
Tp   = 2e-3;
x0   = 4000; %4000;
y0   = 90000;
z0   = 0;
Lam  = 0.03;
dTeta= 0.5;  % 

dev  = 150e6;
dt   = 1/dev;
dxI  = c/(2*dev);
tau  = 6e-6;

Teta = 90 - atan(x0/y0).*gr;
fprintf("Look angle %2.2f \n", Teta);


dl   = 1;
Tsyn = (Lam * sqrt(x0^2 + y0^2 + zsar^2))/(dl*2*Vsar*sin(Teta/gr));
fprintf("Tsyn is %2.2f \n", Tsyn);


My   = 2*ceil(.5*(Tsyn/Tp))
xsar = (-My/2 : My/2-1)*Tp*Vsar;
ty   = (-My/2 : My/2-1)*Tp;

R       = sqrt((x0 - xsar).^2 + y0^2 + zsar^2);
Rcenter = sqrt((x0 - xsar).^2 + y0^2 + zsar^2);

dX = sqrt(x0^2 + y0^2 + zsar^2)*dTeta/gr;
tmin = 2*(sqrt(x0^2 + y0^2 + zsar^2) - dX/2)/c;
tmax = 2*(sqrt(x0^2 + y0^2 + zsar^2) + dX/2)/c + tau;
Mx   = 2*ceil((tmax - tmin)/2/dt)
tx   = tmin + (0:Mx-1)*dt;


s_raw = zeros(My, Mx);
for ny = 1 : My
td          = tx - 2*R(ny)/3e8; 
s_raw(ny, :)= exp(1i*pi*dev/tau*(td.^2-td*tau))*exp(-1i*4*pi*R(ny)./Lam).*(td>=0 & td<=tau); %
end

figure
imagesc(real(s_raw))
title("SAR Raw signal: real part")
grid on


u = s_raw;

%% demodulation
s_demod  = zeros(My, Mx);
Fdd      = 2*Vsar/Lam*cos(Teta/gr);
s_DopMid = exp(-1i*2*pi*Fdd.*ty);
for ny = 1 : My
    s_demod(ny, :) = s_raw(ny, :).*s_DopMid(ny);
end




%% conv in Range direction

R_op    = sqrt(x0^2 + y0^2 + zsar^2);
td0     = tx - 2*R_op/3e8;
h_range = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
hF_range= fft(h_range);
s_range = zeros(My, Mx); 
fs_raw  = zeros(My, Mx); 
fc_raw  = zeros(My, Mx);
for k = 1 : My
     fs_raw(k , :)   = fft(s_demod(k, :));
     fc_raw(k , :)   = fs_raw(k , :).*conj(hF_range);
    s_range(k , :)   = fftshift(ifft(fc_raw(k , :)));
end

figure
imagesc(abs(s_range))
title('Range compression')
xlabel('range time')
ylabel('azimuth time')
grid on



%% range cells correction
fa       =  -1/Tp/2:1/Tsyn:1/Tp/2+1/Tsyn;
% Rcenter;
% dD       =  R_op*(sqrt(fa.^2*Lam^2/(4*Vsar^2)+1)-1);
D  = sqrt(1 - fa.^2*Lam^2/(4*Vsar^2));
dD = R_op.*((1 - D)./D);
% tcen     = Rcenter - R_op;
rangD    =  round(dD/dxI);
rangDmax =  max(rangD);

figure
plot(rangD , '.-b')
ylabel('cell')
grid on

fsmb  = zeros(My, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l))); 
end

fsmbF = fsmb;
smbF  = zeros(My, Mx);


for k=1:My
    for m=1:Mx-rangDmax  
        fsmbF(k,m)=fsmb(k,m+rangD(k));       
    end
end

s_rng = zeros(My, Mx);
for k = 1 : Mx
    s_rng(  :,k)   = fftshift(ifft(fsmbF(:, k )));
end

figure
imagesc(1:Mx, fa, abs(fsmb))
title('Range compression')
xlabel('range time')
ylabel('azimuth frequency')
grid on

figure
imagesc(abs(fsmbF))
title('Range correction')
xlabel('range time')
ylabel('azimuth frequency')
grid on

%% conv in Azimuth direction

R0    = sqrt(y0^2 + zsar^2);
Ka    = 2*Vsar^2/(Lam*R0)*sin(Teta/gr);
smb0  = exp(1i*pi*Ka.*(ty.^2-ty.*Tsyn));     
fsmb0 = fftshift(fft(smb0)); %     
 fsac = zeros(My, Mx); 
  sac = zeros(My, Mx); 
  sF_range = zeros(My, Mx); 
for l = 1 : Mx
sF_range(:,l) =  fftshift(fft(s_range(:,l)));
    fsac(:,l) =  fsmbF(:, l).*(fsmb0).'; %sF_range(:,l)   
     sac(:,l) = (ifft((fsac(:,l)))); %
end

figure
imagesc(abs(sac))
title('Radar Image')
grid on


fprintf("End SAR Model \n");

