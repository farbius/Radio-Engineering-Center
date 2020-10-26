%%
%

clc
clear
close all

fprintf("Start SAR Model \n");
%% SAR parameters
gr   = 180 / pi;
Vsar = 250; 
zsar = 10000;
Tsyn = 5.0;
Tp   = 2e-3;
x0   = 0; %4000;
y0   = 60000;
z0   = 0;
Lam  = 0.03;

dev  = 150e6;
dt   = 1/dev;
tau  = 6e-6;

Teta = 90 - atan(x0/y0).*gr
%% jam parameters
% Type of jamm
% 0 - without any jamm
% 1 - range    0 / pi
% 2 - range    M / pi
% 3 - azimuth  0 / pi
% 4 - azimuth  M / pi

TYPE_Rg = 0; 
TYPE_Az = 0;
T_elr   = 0.06e-6; % elementary pulse duration: range
T_ela   = 10e-3;% elementary pulse duration: azimuth
L_MP    = 8;    % M posl 2^L_MP-1

x_jamm = 0;
y_jamm = 60000;

% [751 1501]

%%
%  look angle

Teta = 90 - atan(x0/y0).*gr;
fprintf("Look angle %2.2f \n", Teta);

My   = 2*ceil(.5*(Tsyn/Tp))
xsar = (-My/2 : My/2-1)*Tp*Vsar;
ty   = (-My/2 : My/2-1)*Tp;

R    = sqrt((x0 + xsar).^2 + y0^2 + zsar^2);

tmin = 400e-6;
tmax = 420e-6;
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

%% jamm signals

% range direction
FAZ_0pi = zeros(1,Mx);

% azimuth direction
FAZ_0piAz = zeros(1,My);
FAZ_MpiAz = zeros(1,My);

jamm      = zeros(My, Mx);
RJ        = sqrt((x_jamm - xsar).^2 + y_jamm^2 + zsar^2);

rect = zeros(1, Mx);

for ny = 1 : My
    
    td = tx - 2*RJ(ny)/3e8;
    rect(td>=0 & td<=tau) = 1; % îãèáàşùàÿ
    jamm(ny, :) = exp(1i*pi*dev/tau*(td.^2-td*tau))*exp(-1i*4*pi*RJ(ny)./Lam).*(td>=0 & td<=tau);
%% TYPE 1    
    if TYPE_Rg == 1
        
        Nelm = ceil (tau / T_elr);
        indx = find(rect == 1);
        FAZ_0pi(indx(1)) = 1;
        IM = 1;
        for k = indx(2) : indx(length(indx))
            if mod(IM,ceil(Mx/Nelm))==0
                    FAZ_0pi(k)=-FAZ_0pi(k-1);
            else
                    FAZ_0pi(k)= FAZ_0pi(k-1);
            end
            IM = IM + 1;
        end
        jamm(ny, :) = jamm(ny, :).*FAZ_0pi.*(td>=0 & td<=tau);
    end % TYPE 1
%% TYPE 2
       if TYPE_Rg == 2
            NM           = 2^L_MP-1;   
            MP           = M_POSL(L_MP,0);
            indx         = find(rect == 1);
            iM           = 1; 
            IM = 1;
            FAZ_0pi(indx(1)) = MP(iM);
            
            for k = indx(2) : indx(length(indx))
                if mod(IM,ceil(length(indx)/NM))==0
                        if iM<NM
                          iM=iM+1;
                        end
                        FAZ_0pi(k)=MP(iM);
                else
                        FAZ_0pi(k)= FAZ_0pi(k-1);
                end
                IM = IM + 1;
            end
            jamm(ny, :) = jamm(ny, :).*FAZ_0pi.*(td>=0 & td<=tau);
        end % TYPE 2
%% TYPE 3
        if TYPE_Az == 1 %0/pi 
            Ni=floor(Tsint/T_ela); %    
            Ki=floor(My/Ni); % 
            FAZ_0piAz(1)=1;
            for ki=2:My
                if mod(ki,Ki)==0
                    FAZ_0piAz(ki)=-FAZ_0piAz(ki-1);
                else
                    FAZ_0piAz(ki)= FAZ_0piAz(ki-1);
                end
            end
            jamm(ny,:)=jamm(ny,:)*FAZ_0piA(ny).*(td>=0 & td<=tau);
        end % TYPE 3
%% TYPE 4
        if TYPE_Az == 2 %
            Ni=floor(Tsint/T0piJA); %     
            Ki=floor(My/Ni); %
            iM=1;  FAZ_MpiAz(1)=MP(iM);
            for ki=2:My
                if mod(ki,Ki)==0
                    if iM<NM
                        iM=iM+1;
                    end
                    FAZ_MpiAz(ki)=MP(iM);
 
                else
                    FAZ_MpiAz(ki)=FAZ_MpiAz(ki-1);
                end
            end
            jamm(ny,:)=jamm(ny,:)*FAZ_MpiAz(ny).*(td>=0 & td<=tau);
        end %TypeJA=3
  
end

u = s_raw + 5.*jamm;




%% conv in Range direction

R_op    = sqrt(x0^2 + y0^2 + zsar^2);
td0     = tx - 2*R_op/3e8;
h_range = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
hF_range= fft(h_range);
s_range = zeros(My, Mx); sj_range = zeros(My, Mx);
fs_raw  = zeros(My, Mx); fsj_raw  = zeros(My, Mx);
fc_raw  = zeros(My, Mx); fcj_raw  = zeros(My, Mx);
for k = 1 : My
     fs_raw(k , :)   = fft(1.5.*jamm(k, :));
     fsj_raw(k , :)  = fft(s_raw(k, :));
     fc_raw(k , :)   = fs_raw(k , :).*conj(hF_range);
     fcj_raw(k , :)  = fsj_raw(k , :).*conj(hF_range);
    s_range(k , :)   = fftshift(ifft(fc_raw(k , :)));
   sj_range(k , :)   = fftshift(ifft(fcj_raw(k , :)));
end

figure
imagesc(abs(s_range))
title('Range compression')
grid on

figure
plot(tx./1e-6, abs(s_range(751, :)), '.-b')
xlim([405 415])
t = title('Âûõîä ñîãëàñîâàííîãî ôèëüòğà: êàíàë äàëüíîñòè');
t.FontSize = 13;
x = xlabel('Âğåìÿ, ìêñ');
x.FontSize = 12;
grid on
%% pictures
figure
subplot(2,1,1)
plot(tx./1e-6, real(jamm(751, :)), '.-b', tx./1e-6, abs(jamm(751, :)), '.-k')
ylim([-1.5 1.5])
xlim([405 412])
t = title('Ôàçîìàíèïóëèğîâàííûé èìïóëüñ');
t.FontSize = 13;
x = xlabel('Âğåìÿ, ìêñ');
x.FontSize = 12;
grid on
subplot(2,1,2)
plot(tx./1e-6, FAZ_0pi, '.-b')
ylim([-1.5 1.5])
xlim([405 412])
t = title('Çàêîí ôàçîâîé ìàíèïóëÿöèè');
t.FontSize = 13;
x = xlabel('Âğåìÿ, ìêñ');
x.FontSize = 12;
grid on

%% conv in Azimuth direction

R0    = sqrt(y0^2 + zsar^2);
Ka    = 2*Vsar^2/(Lam*R0);
smb0  = exp(1i*pi*Ka.*(ty.^2-ty.*Tsyn));     
fsmb0 = fftshift(fft(smb0)); %     
 fsac = zeros(My, Mx); 
  sac = zeros(My, Mx); 
  sF_range = zeros(My, Mx); 
for l = 1 : Mx
sF_range(:,l) =  fftshift(fft(s_range(:,l)));
    fsac(:,l) =  sF_range(:,l).*(fsmb0).'; %   
     sac(:,l) = (ifft((fsac(:,l)))); %
end

figure
imagesc(abs(sac))
title('Radar image')
grid on


fprintf("End SAR Model \n");

