%%
%
clc
clear
close all
LOAD = 1;
fprintf(">> Start SAR Model \n");
%% константы
gr   = 180 / pi;
c    = 3e8;
%% параметры БРЛС
Vsar = 250;   % скорость движения БРЛС
zsar = 10000; % высота полета
% центр участка картографирования
x0   = 10000; 
y0   = 0;
z0   = 0;

fc    = 10e9;        % частота зондирующего сигнала
Lam   = c/fc;        % длина волны
La    = 1.0;         % ширина раскрыва апертуры реальной ДНА
Teta05= Lam/La*gr;   % ширина ДНА по уровню 0.5
fprintf(">> ширина гл луча ДНА %2.2f град \n", Teta05);

TetaQ = 5; %
fprintf(">> угол наклона гл луча ДНА %2.2f град\n", TetaQ);

%% зондирующий сигнал
dev  = 150e6;  % девиация ЛЧМ  импульса
tau  = 6e-6;   % длительность ЛЧМ импульса
dt   = 1/dev/2;% период работы АЦП
dxI  = c/dev;  % шаг по дальности
dyI  = La/2;   % шаг по азимуту



%% формирование матрицы целевой обстановки
targets = rgb2gray(imread('targets', 'png'));
[Nx, Ny] = size(targets);
Ntarget = Nx*Ny;
num=1; xn=zeros(Ntarget,1); 
yn=xn; Fn=xn; 
% формирование векторов координат участка картографирования
% координата z отсутствует 
for iNy=1:Ny 
    for iNx=1:Nx 
        xn(num)=(iNx-Nx/2);   
        yn(num)=(Ny/2-iNy+1);  
        Fn(num)=double(targets(iNy,iNx))/1;     
        num=num+1; 
    end 
end 
% коэффициент растяжения
kx = 3; ky = 3;
xn=xn*kx*dxI; yn=yn*ky*dyI;


%% наклонная дальность центра участка картографирования и положение БРЛС
R_0  = sqrt(x0^2 + zsar^2);
% положение БРЛС на интервале синтезирования
u1   = -R_0*tan((TetaQ+.5*Teta05)/gr); % начало интервала
u2   = -R_0*tan((TetaQ)/gr);           % середина интервала
u3   = -R_0*tan((TetaQ-.5*Teta05)/gr); % окончание интервала
fprintf(">> положение БРЛС по оси u\n");
fprintf(">> u1 = %5.2f м, u2 = %5.2f м, u3 = %5.2f м,\n", u1, u2, u3);
% наклонные дальности на интервале синтезтрования
R3   = sqrt(u3.^2   + R_0^2);
R2   = sqrt(u2.^2   + R_0^2);
R1   = sqrt(u1.^2   + R_0^2);
fprintf(">> наклонная дальность R1 = %8.2f m, R2 = %8.2f m, R3 = %8.2f m\n", R1, R2, R3);
% Доплер на интервале синтезирования
fD1  = 2*Vsar/Lam*sin((TetaQ-.5*Teta05)/gr);
fD3  = 2*Vsar/Lam*sin((TetaQ+.5*Teta05)/gr);
fD   = fD3 - fD1;
fD2  = fD1 + fD/2;
fprintf(">> нижний Доплер fD1 = %8.2f Гц \n", fD1);
fprintf(">> верхний Доплер fD3 = %8.2f Гц \n", fD3);
fprintf(">> средний Доплер fD2 = %8.2f Гц \n", fD2);
fprintf(">> ширина спектра тр сигнала %8.2f Гц \n", fD);

%% TODO
% рассчитать период повторения
% для переноса спектра траекторного сигнала 
% в первую зону Найквиста
Fprf = 700;
Tp   = 1/Fprf;

if (Fprf > c/(2*abs(R1)))
    fprintf(">> неоднозн по дальности \n");
    return
elseif(Fprf < 2*Vsar/La)
    fprintf(">> неоднозн по Доплеру \n");  
    return
end

fprintf(">> выбор частоты повторения Fprf\n");
fprintf(">> %8.2f <= %d <= %8.2f \n", 2*Vsar/La, Fprf, c/(2*abs(R1)));

%% Формирование разверток по азимуту и дальности
% ширина участка синтезирования
Ls   = R_0*(tan((TetaQ+.5*Teta05)/gr) - tan((TetaQ-.5*Teta05)/gr));% length aperture
fprintf(">> ширина участка синтезирования %5.2f м \n", Ls);
% длительность интервала синтезирования
Tsyn = Ls/Vsar;
fprintf(">> длительность интервала синтезирования %2.2f с \n", Tsyn);
% колличество накопленных импульсов
My   =   ceil(Tsyn/Tp);
% развертка "медленного" времени
ty   = (0 : My-1)*Tp;
% время начала и окончания приема отраженного сигнала
tmin = 2*R_0/c - 0.5*tau;
tmax = 2*R_0/c + 1.5*tau;
% количество отсчетов по дальности
Mx   = 2*ceil((tmax - tmin)/2/dt);
% развертка "быстрого" времени
tx   = tmin + (0:Mx-1)*dt;
fprintf(">> матрица РСА Mx= %d, My= %d \n", Mx, My);
%%
% изменение наклонной дальности
u = linspace(u1, u3, My);
R_center   = sqrt(u.^2 + R_0^2);
Teta = linspace(TetaQ+.5*Teta05, TetaQ-.5*Teta05, My);
f_dop= 2*Vsar/Lam*sin(Teta./gr);

figure
plot(R_center, f_dop, '.-b', R3, fD1, 'or', R2, fD2, 'or', R1, fD3, 'or')
txt = ['Teta_H: ' num2str(TetaQ) ' град'];
text(R2+1.0,fD2,txt)
title('Изменение Доплера на интервале синтезирования')
xlabel('Range, m')
ylabel('Доплер, Гц')
grid on

%% формирование отраженного сигнала (матрицы РСА)
% координата u БРЛС на интервале синтезирования
u = linspace(u1, u3, My);
% матрица РСА
s_raw = zeros(My, Mx);
if LOAD == 1
    load('s_raw.mat');
else
for ny = 1 : My
    for k = 1 : Ntarget
          R = sqrt((u(ny) + yn(k)).^2 + (x0 + xn(k))^2 + (zsar + z0)^2);
         td = tx - 2*R/3e8; 
s_raw(ny, :)= s_raw(ny, :) + Fn(k).*exp(1i*pi*dev/tau*(td.^2-td*tau))*exp(-1i*4*pi*R./Lam).*(td>=0 & td<=tau);
    end
    if(mod(ny, 100)==0)
    fprintf("...%d", round(ny/My*100));
    end
end

fprintf("/n");
end

fprintf("\n");

figure
imagesc(real(s_raw))
title("SAR Raw signal: real part")
grid on

%% conv in Range direction
td0      = tx - 2*R_0/3e8;
h_range  = exp(1i*pi*dev/tau*(td0.^2-td0*tau)).*(td0>=0 & td0<=tau);
hF_range = fft(h_range);
s_range  = zeros(My, Mx); 
fs_raw   = zeros(My, Mx);
fc_raw   = zeros(My, Mx);
for k = 1 : My
     fs_raw(k , :)   = fft(s_raw(k, :));
     fc_raw(k , :)   = fs_raw(k , :).*conj(hF_range);
    s_range(k , :)   = fftshift(ifft(fc_raw(k , :)));
end


figure
imagesc(1:Mx, 1:My, abs(s_range))
title('Range compression')
xlabel('range time bins')
ylabel('azimuth time bins')
grid on

%% range cells correction
% FFT length
NAzFFT = 1024;
fa     = (1 : NAzFFT);
% БПФ по азимуту
fsmb  = zeros(NAzFFT, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l), NAzFFT)); 
end
% Функция для коррекции 
N_1 = R3*Lam^2/(4*Vsar^2/dxI);
dNr = ceil(N_1.*(fD1 - fa./(Tp*NAzFFT)).^2 - N_1*fD1^2);
% broadside view
% dNr = R_0*sqrt(1 - (k.*Lam*Fprf/(2*Vsar*NAzFFT)).^2) - R_0;
% коррекция миграции дальности
fsmbF = fsmb;
smbF  = zeros(NAzFFT, Mx);
for k=1:NAzFFT  
        fsmbF(k,:)=circshift(fsmb(k,:), dNr(k));       
end

figure
plot( fa,dNr, '.-b')
xlabel('range bins')
ylabel('azimuth freq bins')
title('Range migration')
grid on

figure
imagesc(1:Mx, fa, abs(fsmb))
title('Range compression: Az FFT ')
xlabel('range time bins')
ylabel('azimuth freq bins')
grid on

figure
imagesc(1:Mx, fa, abs(fsmbF))
title('Range correction: Az FFT')
xlabel('range time bins')
ylabel('azimuth freq bins')
grid on


%% conv in Azimuth direction
% крутизна траекторного сигнала
Ka    = 2*Vsar^2/(Lam*R_0);
% фаза траекторного сигнала
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
ylabel('Dopler time bins')
grid on

figure
mesh(abs(sac))
grid on

fprintf(">> End SAR Model \n");

