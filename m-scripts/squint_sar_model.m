%%
%

clc
clear
close all

fprintf("Start SAR Model \n");
%% константы
gr   = 180 / pi;
c    = 3e8;
%% параметры БРЛС
Vsar = 250; 
zsar = 10000;

% центр участка картографирования
x0   = 10000; 
y0   = 800;
z0   = 0;

x1   = 10000;
y1   = 810;
z1   = 0;

fc    = 10e9;        % частота зондирующего сигнала
Lam   = c/fc;        % длина волны
La    = 1.0;         % ширина раскрыва апертуры реальной ДНА
Teta05= Lam/La*gr;   % ширина ДНА по уровню 0.5
fprintf(">> ширина гл луча ДНА %2.2f град \n", Teta05);

TetaQ = atan(y0/x0).*gr;
fprintf(">> угол наклона гл луча ДНА %2.2f град\n", TetaQ);

%% зондирующий сигнал
dev  = 150e6;
dt   = 1/dev/2;
dxI  = c/dev;
dl   = 1;
tau  = 6e-6;

%% TODO
% рассчитать период повторения
% для переноса спектра тр сигнала 
% в первую зону Найквиста

Fprf = 700;
Tp   = 1/Fprf;
%%

Ls   = sqrt(x0^2 + zsar^2)*(tan((TetaQ+.5*Teta05)/gr) - tan((TetaQ-.5*Teta05)/gr));% length aperture
fprintf(">> ширина участка синтезирования %5.2f м \n", Ls);
R_0  = sqrt(0^2 + x0^2 + zsar^2);
%% координата x БРЛС в начале, середине и конце интервала синтезирования
Tsyn = Ls/Vsar;
fprintf(">> длительность интервала синтезирования %2.2f с \n", Tsyn);
My   =   ceil(Tsyn/Tp);
ty   = (0 : My-1)*Tp;
fa   = (0 : My-1)/Tsyn;


u1   = -R_0*tan((TetaQ+.5*Teta05)/gr); % 
u2   = -R_0*tan((TetaQ)/gr); %
u3   = -R_0*tan((TetaQ-.5*Teta05)/gr); % 
fprintf(">> положение БРЛС по оси X\n");
fprintf(">> x1 = %5.2f м, x2 = %5.2f м, x3 = %5.2f м,\n", u1, u2, u3);


fD1  = 2*Vsar/Lam*sin((TetaQ-.5*Teta05)/gr);
fprintf(">> нижний Доплер fD1 = %8.2f Гц \n", fD1);
fD3  = 2*Vsar/Lam*sin((TetaQ+.5*Teta05)/gr);
fprintf(">> верхний Доплер fD3 = %8.2f Гц \n", fD3);
fD   = fD3 - fD1;
fD2  = fD1 + fD/2;
fprintf(">> средний Доплер fD2 = %8.2f Гц \n", fD2);
fprintf(">> ширина спектра тр сигнала %8.2f Гц \n", fD);


R3   = sqrt(u3.^2   + R_0^2);
R2   = sqrt(u2.^2   + R_0^2);
R1   = sqrt(u1.^2   + R_0^2);
fprintf(">> slant range R1 = %8.2f m, R2 = %8.2f m, R3 = %8.2f m\n", R1, R2, R3);

if (Fprf > c/(2*abs(R1)))
    fprintf(">> неоднозн по дальности \n");
    return
elseif(Fprf < 2*Vsar/La)
    fprintf(">> неоднозн по Доплеру \n");  
    return
end

fprintf(">> выбор частоты повторения Fprf\n");
fprintf(">> %8.2f <= %d <= %8.2f \n", 2*Vsar/La, Fprf, c/(2*abs(R1)));

%% TODO 
% пересчитать развертку по дальности 
% из кинематики БРЛС
dX   =   sqrt(x0^2 + y0^2 + zsar^2)*Teta05/gr;
tmin = 2*R_0/c - 0.5*tau;
tmax = 2*R_0/c + 1.5*tau;
Mx   = 2*ceil((tmax - tmin)/2/dt);
 

fprintf(">> матрица РСА Mx= %d, My= %d \n", Mx, My);

tx   = tmin + (0:Mx-1)*dt;

% координата X БРЛС
Xsar = linspace(u1, u3, My);
% изменение наклонной дальности
R    = sqrt(Xsar.^2 + R_0^2);


Teta = linspace(TetaQ+.5*Teta05, TetaQ-.5*Teta05, My);
f_dop= 2*Vsar/Lam*sin(Teta./gr);

figure
plot(R, f_dop, '.-b', R3, fD1, 'or', R2, fD2, 'or', R1, fD3, 'or')
txt = ['Teta_H: ' num2str(TetaQ) ' град'];
text(R2+1.0,fD2,txt)
title('Изменение Доплера на интервале синтезирования')
xlabel('Дальность, м')
ylabel('Доплер, Гц')
grid on


s_raw = zeros(My, Mx);
for ny = 1 : My
td          = tx - 2*R(ny)/3e8; 
s_raw(ny, :)= exp(1i*pi*dev/tau*(td.^2-td*tau))*exp(-1i*4*pi*R(ny)./Lam).*(td>=0 & td<=tau);
end


figure
imagesc(real(s_raw))
title("SAR Raw signal: real part")
grid on


%% demodulation
s_demod  = zeros(My, Mx);
s_DopMid = exp(1i*2*pi*fD2.*ty);
for ny = 1 : My
    s_demod(ny, :) = s_raw(ny, :); %.*s_DopMid(ny);
end



%% conv in Range direction

td0      = tx - 2*R_0/3e8;
h_range  = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
hF_range = fft(h_range);
s_range  = zeros(My, Mx); fs_raw  = zeros(My, Mx);fc_raw  = zeros(My, Mx);
for k = 1 : My
     fs_raw(k , :)   = fft(s_demod(k, :));
     fc_raw(k , :)   = fs_raw(k , :).*conj(hF_range);
    s_range(k , :)   = fftshift(ifft(fc_raw(k , :)));
end
% 
figure
imagesc(1:Mx, 1:My, abs(s_range))
title('Range compression')
xlabel('range time bins')
ylabel('azimuth time bins')
grid on



%% range cells correction
% FFT length
NAzFFT = 1024;
%% TODO
% пересчитать спектр траекторного сигнала
% с учетом субдискретизации
fa = (1 : NAzFFT);

% cells migration
dN  = R1 - R3;

% squint view
N_1 = R3*Lam^2/(4*Vsar^2/dxI);
dNr = ceil(N_1.*(fD1 - fa./(Tp*NAzFFT)).^2 - N_1*fD1^2);

% broadside view
% dNr = R_0*sqrt(1 - (k.*Lam*Fprf/(2*Vsar*NAzFFT)).^2) - R_0;

figure
plot( fa,dNr, '.-b')
xlabel('range bins')
ylabel('azimuth freq bins')
title('Range migration')
grid on

fsmb  = zeros(NAzFFT, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l), NAzFFT)); 
end


figure
imagesc(1:Mx, fa, abs(fsmb))
title('Range compression: Az FFT ')
xlabel('range time bins')
ylabel('azimuth freq bins')
grid on

fsmbF = fsmb;
smbF  = zeros(NAzFFT, Mx);
for k=1:NAzFFT  
        fsmbF(k,:)=circshift(fsmb(k,:), dNr(k));       
end



figure
% imagesc(tx.*c/2./1e3, fa, abs(fsmbF))
imagesc(1:Mx, fa, abs(fsmbF))
title('Range correction: Az FFT')
xlabel('range time bins')
ylabel('azimuth freq bins')
grid on


%% conv in Azimuth direction

Ka    = 2*Vsar^2/(Lam*R_0);
% smb0  = exp(1i*pi*Ka.*(ty.^2-ty.*Tsyn));
FazOp = 1i*2*pi*fD2.*ty + 1i*pi*Ka.*ty.^2;
smb0  = exp(FazOp);
fsmb0 = fftshift(fft(smb0, NAzFFT)); %     
 fsac = zeros(NAzFFT, Mx); 
  sac = zeros(NAzFFT, Mx); 
  sF_range = zeros(NAzFFT, Mx); 
for l = 1 : Mx
    fsac(:,l) =  fsmbF(:, l).*fsmb0.';  
     sac(:,l) =  ifft(fsac(:,l)); 
end

figure
imagesc(1:Mx, fa, abs(sac))
title('Radar Image')
xlabel('range time bins')
ylabel('Dopler freq bins')
grid on

figure
mesh(abs(sac))
grid on





% My = 4096;
% 
% fa       =  (-My/2 : My/2-1)/Tsyn;
% fa = fa + 2*Vsar/Lam*cos((Teta)/gr); %fa + 2*Vsar/Lam*cos((Teta)/gr)
% % dD_w       =  R_op*(sqrt(fa.^2*Lam^2/(4*Vsar^2)+1)-1);
% dD_w = R_op.*sqrt(1 - (Lam^2.*fa.^2)./(4*Vsar^2)) - R_op;
% rangD    =  ceil(dD_w/dxI);
% 
% fDL      = 2*Vsar/Lam*cos((Teta - Teta05/2)/gr) - 2*Vsar/Lam*cos((Teta)/gr);
% rangDw   = R_op*Lam^2/(4*Vsar^2*dxI).*(fDL - fa) - R_op*Lam^2*fDL^2/(4*Vsar^2*dxI);
% rangDw   = ceil(rangDw);
% 
% figure
% plot(rangD , '.-b')
% ylabel('cell')
% grid on
% 
% fsmb  = zeros(My, Mx);
% parfor l=1:Mx
%     fsmb(:,l)=fftshift(fft(s_range(:,l), My)); 
% end
% 
% figure
% imagesc(abs(fsmb))
% title('Range compression')
% xlabel('range time')
% ylabel('azimuth frequency')
% grid on
% 
% 
% fsmbF = fsmb;
% smbF  = zeros(My, Mx);
% 

% 

% 


fprintf("End SAR Model \n");

