%% алгоритм формирования принятого сигнала 
% Ростов Алексей, 17.04.20
% farbius@protonmail.com
clc
clear
close all
fprintf('> Start script \n\r');
%%  I этап
fprintf('> Загрузка предустановленных параметров ...\n\r');
init_parameters;
fprintf('> Done \n\r');
%% Формирование разверток матричного сигнала РСА
  
fprintf('> Расчет векторов slow и fast time...\n\r');
calc_timevector;
fprintf('> Done \n\r');

%% Формирование подстилающей поверхности (ПП) %%%%%%%%%
% не входит в алгоритм формирования принятого сигнала %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('> Формирование подстилающей поверхности ...\n\r');
calc_terrain;
fprintf('> Done \n\r');

%% Формирование матрицы внутренних шумов %%%%%%%%%%%%%%              
% не входит в алгоритм формирования принятого сигнала %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('> Формирование матрицы внутренних шумов ...\n\r');
calc_noise;
fprintf('> Done \n\r'); 

%% Формирование матрицы активных помех %%%%%%%%%%%%%%%%              
% не входит в алгоритм формирования принятого сигнала %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('> Формирование матрицы активных помех ...\n\r');
calc_jam;
fprintf('> Done \n\r'); 

%% II этап
fprintf('> Расчет матрицы РСА ...\n\r');
s = zeros(My, Mx); % матрица сигнала RAW сигнала SAR
LOAD_MAT = 1;

if LOAD_MAT == 1
   load('s_raw.mat'); 
else
    parfor ny = 1 : My
        for m = 1 : Num_trg
           % расчет дальности точки в нулевой момент времени
           r0    = sqrt(x_p(m)^2 + y_p(m)^2 + z0^2);
           % расчет текущей дальности
           R     = sqrt(r0^2 - 2*x_p(m)*V*t_y(ny) + (V*t_y(ny))^2); 
           % расчет времени задержки сигнала
%            td    = t_x - 2*R_0/c;
           td    = t_x - 2*R_0/c;
           s(ny,:) = s(ny,:)+Sigma(m).*exp(-1i*(4*pi*R/Lam-pi*dev/tau*(td.^2-td*tau))).*(td>=0 & td<=tau);
        end
    end
    save('s_raw.mat','s'); 
end % LOAD MAT
fprintf('> Done \n\r');

%% III этап
% аддитивный сигнал (сигнал, помеха, шум)
u = s;
% u = s + jam + n;

figure(1)
imagesc( t_x./1e-6, t_y, real(u))
colormap jet
title('Матрица РСА: реальная квадратура')
xlabel('fast time, us')
ylabel('slow time, s') 
grid on

figure(2)
imagesc( t_x./1e-6, t_y, imag(u))
colormap jet
title('Матрица РСА: мнимая квадратура')
xlabel('fast time, us')
ylabel('slow time, s') 
grid on

figure(3)
imagesc( t_x./1e-6, t_y, abs(u))
colormap jet
title('Матрица РСА: огибающая')
xlabel('fast time, us')
ylabel('slow time, s') 
grid on

figure(4)
plot( t_x./1e-6, abs(u(1, :)))
title('Принятый импульс РСА: огибающая')
xlabel('fast time, us')
grid on

figure(5)
plot( t_x./1e-6, real(u(1, :)))
title('Принятый импульс РСА: real')
xlabel('fast time, us')
grid on

figure(6)
plot( t_x./1e-6, imag(u(1, :)))
title('Принятый импульс РСА: imag')
xlabel('fast time, us')
grid on
% 
% figure(7)
% plot( t_y, imag(u(:, 4500)))
% title('Траекторный сигнал РСА: imag')
% xlabel('fast time, us')
% grid on
  