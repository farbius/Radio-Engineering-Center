%%
%

clc
clear
close all

gr   = 180 / pi;
c    = 3e8;
Vsar = 250; 
zsar = 10000;
Tp   = 2e-3;
x0   = 4000;
y0   = 90000;
Lam  = 0.03;

dev  = 150e6;
dt   = 1/dev;
dxI  = c/(2*dev);
tau  = 6e-6;
dTeta= 0.5;  

dX = sqrt(x0^2 + y0^2 + zsar^2)*dTeta/gr;
tmin = 2*(sqrt(x0^2 + y0^2 + zsar^2) - dX/2)/c;
tmax = 2*(sqrt(x0^2 + y0^2 + zsar^2) + dX/2)/c + tau;
Mx   = 2*ceil((tmax - tmin)/2/dt)
tx   = tmin + (0:Mx-1)*dt;


Teta = 90 - atan(x0/y0).*gr;
fprintf("Look angle %2.2f \n", Teta);

dl    = 1;
Tsyn    = (Lam * sqrt(x0^2 + y0^2 + zsar^2))/(dl*2*Vsar*sin(Teta/gr));
My      = 2*ceil(.5*(Tsyn/Tp))
ty      = (-My/2 : My/2-1)*Tp;
fa       =  -1/Tp/2:1/Tsyn:1/Tp/2+1/Tsyn;


R0      = sqrt(x0.^2 + y0^2 + zsar^2);
R       = sqrt((x0 - Vsar.*ty).^2 + y0^2 + zsar^2);

f_n     =(R0*Lam)/(2*Vsar^2*sin(Teta/gr)).*fa;
R_f     = sqrt((x0 - Vsar.*f_n).^2 + y0^2 + zsar^2);
load('s_range.mat');
D    = sqrt(1 - fa.^2*Lam^2/(4*Vsar^2));
D_sq = R_f./R0;

dD    = R0.*((1 - D)./D);
dD_sq = R0.*((1 - D_sq)./D_sq);

rangD    =  round(dD/dxI);
randD_sq =  round(dD_sq/dxI);

rangDmax =  max(abs(randD_sq));

figure
plot(randD_sq , '.-b')
ylabel('cell')
xlabel('range cells')
grid on

fsmb  = zeros(My, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l))); 
end

fsmbF = fsmb;
smbF  = zeros(My, Mx);

figure
imagesc(1:Mx, fa, abs(fsmb))
title('Range compression')
xlabel('range cells')
ylabel('azimuth frequency')
grid on



for k=1:My
%     for m=1:Mx  
        fsmbF(k,:) = circshift(fsmb(k, :), -1*randD_sq(k));
%         if(randD_sq(k) >= 0)
%         fsmbF(k,m)=fsmb(k,m+randD_sq(k));  
%         else
%         fsmbF(k,m)=fsmb(k,m-randD_sq(k));    
%         end
%     end
end

s_rng = zeros(My, Mx);
for k = 1 : Mx
    s_rng(  :,k)   = fftshift(ifft(fsmbF(:, k )));
end


figure
imagesc(1:Mx, fa, abs(fsmbF))
title('Range correction')
xlabel('range cells')
ylabel('azimuth frequency')
grid on

figure
imagesc(1:Mx, fa, abs(s_rng))
title('Range correction')
xlabel('range cells')
ylabel('azimuth frequency')
grid on


% 
% R_t     = R0 + (Vsar.*ty).^2./(2*R0);
% 
% Ka      = 2*Vsar^2/R0/Lam;
% f_a     = -Ka.*ty;
% 
% R_f     = R0 + (R0*Lam^2)/(8*Vsar^2).*f_a.^2;
% 
% dR      = round(R_f./dr);
% 
% figure
% plot(R_t,'.-b')
% hold on
% plot(R_f, 'o-r')
% grid on
% 
% figure
% plot(dR,'.-b')
% grid on