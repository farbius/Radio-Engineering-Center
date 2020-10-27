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
x0   = 1000; %4000;
y0   = 90000;
z0   = 0;
Lam  = 0.03;
Teta05= 0.8;  % 

dev  = 150e6;
dt   = 1/dev/2;
dxI  = c/(2*dev)/2;
tau  = 6e-6;

Teta = 90 - atan(x0/y0).*gr;
fprintf("Look angle %2.2f \n", Teta);


dl   = 1;
Tsyn = (Lam * sqrt(x0^2 + y0^2 + zsar^2))/(dl*2*Vsar*sin(Teta/gr));
fprintf("Tsyn is %2.2f \n", Tsyn);
R_op    = sqrt(x0^2 + y0^2 + zsar^2)
Tp =  1e-3;
fprintf("Pulse repetition %2.2f ms\n", Tp/1e-3);



My   = 2*ceil(.5*(Tsyn/Tp))
% ty   = (-My/2 : My/2-1)*Tp;
ty   = (0:My-1)*Tp;

R       = sqrt(((ty-Tsyn/2 - x0/Vsar).*Vsar).^2 + y0^2 + zsar^2);

dX = sqrt(x0^2 + y0^2 + zsar^2)*Teta05/gr;
tmin = 2*(sqrt(x0^2 + y0^2 + zsar^2) - dX/2)/c;
tmax = 2*(sqrt(x0^2 + y0^2 + zsar^2) + dX/2)/c + tau;
Mx   = 2*ceil((tmax - tmin)/2/dt)
tx   = tmin + (0:Mx-1)*dt;


s_raw = zeros(My, Mx);
parfor ny = 1 : My
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
parfor ny = 1 : My
    s_demod(ny, :) = s_raw(ny, :).*s_DopMid(ny);
end



%% conv in Range direction


td0     = tx - 2*R_op/3e8;
h_range = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
hF_range= fft(h_range);
s_range = zeros(My, Mx); 
fs_raw  = zeros(My, Mx); 
fc_raw  = zeros(My, Mx);
parfor k = 1 : My
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
% My = 4096;

fa       =  (-My/2 : My/2-1)/Tsyn;
fa = fa + 2*Vsar/Lam*cos((Teta)/gr); %fa + 2*Vsar/Lam*cos((Teta)/gr)
% dD_w       =  R_op*(sqrt(fa.^2*Lam^2/(4*Vsar^2)+1)-1);
dD_w = R_op.*sqrt(1 - (Lam^2.*fa.^2)./(4*Vsar^2)) - R_op;
rangD    =  ceil(dD_w/dxI);

fDL      = 2*Vsar/Lam*cos((Teta - Teta05/2)/gr) - 2*Vsar/Lam*cos((Teta)/gr);
rangDw   = R_op*Lam^2/(4*Vsar^2*dxI).*(fDL - fa) - R_op*Lam^2*fDL^2/(4*Vsar^2*dxI);
rangDw   = ceil(rangDw);

figure
plot(rangD , '.-b')
ylabel('cell')
grid on

fsmb  = zeros(My, Mx);
parfor l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l), My)); 
end

figure
imagesc(abs(fsmb))
title('Range compression')
xlabel('range time')
ylabel('azimuth frequency')
grid on


fsmbF = fsmb;
smbF  = zeros(My, Mx);

for k=1:My  
%     for m=350:Mx-350  
        fsmbF(k,:)=circshift(fsmb(k,:), rangD(k));%fsmb(k,m+rangD(k));       
%     end
end

s_rng = zeros(My, Mx);
for k = 1 : Mx
    s_rng(  :,k)   = fftshift(ifft(fsmbF(:, k )));
end



figure
imagesc(abs(fsmbF))
title('Range correction')
xlabel('range time')
ylabel('azimuth frequency')
grid on

%% conv in Azimuth direction
ty   = (-My/2 : My/2-1)*Tp;
R0    = sqrt(y0^2 + zsar^2);
Ka    = 2*Vsar^2/(Lam*R0)*sin(Teta/gr);
smb0  = exp(1i*pi*Ka.*(ty.^2-ty.*Tsyn));     
fsmb0 = fftshift(fft(smb0)); %     
 fsac = zeros(My, Mx); 
  sac = zeros(My, Mx); 
  sF_range = zeros(My, Mx); 
for l = 1 : Mx
    fsac(:,l) =  fsmbF(:, l).*(fsmb0).'; %sF_range(:,l)   
     sac(:,l) = (ifft((fsac(:,l)))); %
end

figure
imagesc(tx.*c/2, ty, abs(sac))
title('Radar Image')
grid on

figure
mesh(abs(sac))
grid on


fprintf("End SAR Model \n");

