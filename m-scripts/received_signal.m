%% �������� ������������ ��������� ������� 
% ������ �������, 17.04.20
% farbius@protonmail.com
clc
clear
close all
fprintf('> Start script \n\r');
%%  I ����
fprintf('> �������� ����������������� ���������� ...\n\r');
init_parameters;
fprintf('> Done \n\r');
%% ������������ ��������� ���������� ������� ���
  
fprintf('> ������ �������� slow � fast time...\n\r');
calc_timevector;
fprintf('> Done \n\r');

%% ������������ ������������ ����������� (��) %%%%%%%%%
% �� ������ � �������� ������������ ��������� ������� %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('> ������������ ������������ ����������� ...\n\r');
calc_terrain;
fprintf('> Done \n\r');

%% ������������ ������� ���������� ����� %%%%%%%%%%%%%%              
% �� ������ � �������� ������������ ��������� ������� %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('> ������������ ������� ���������� ����� ...\n\r');
calc_noise;
fprintf('> Done \n\r'); 

%% ������������ ������� �������� ����� %%%%%%%%%%%%%%%%              
% �� ������ � �������� ������������ ��������� ������� %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('> ������������ ������� �������� ����� ...\n\r');
calc_jam;
fprintf('> Done \n\r'); 

%% II ����
fprintf('> ������ ������� ��� ...\n\r');
s = zeros(My, Mx); % ������� ������� RAW ������� SAR
LOAD_MAT = 1;

if LOAD_MAT == 1
   load('s_raw.mat'); 
else
    parfor ny = 1 : My
        for m = 1 : Num_trg
           % ������ ��������� ����� � ������� ������ �������
           r0    = sqrt(x_p(m)^2 + y_p(m)^2 + z0^2);
           % ������ ������� ���������
           R     = sqrt(r0^2 - 2*x_p(m)*V*t_y(ny) + (V*t_y(ny))^2); 
           % ������ ������� �������� �������
%            td    = t_x - 2*R_0/c;
           td    = t_x - 2*R_0/c;
           s(ny,:) = s(ny,:)+Sigma(m).*exp(-1i*(4*pi*R/Lam-pi*dev/tau*(td.^2-td*tau))).*(td>=0 & td<=tau);
        end
    end
    save('s_raw.mat','s'); 
end % LOAD MAT
fprintf('> Done \n\r');

%% III ����
% ���������� ������ (������, ������, ���)
u = s;
% u = s + jam + n;

figure(1)
imagesc( t_x./1e-6, t_y, real(u))
colormap jet
title('������� ���: �������� ����������')
xlabel('fast time, us')
ylabel('slow time, s') 
grid on

figure(2)
imagesc( t_x./1e-6, t_y, imag(u))
colormap jet
title('������� ���: ������ ����������')
xlabel('fast time, us')
ylabel('slow time, s') 
grid on

figure(3)
imagesc( t_x./1e-6, t_y, abs(u))
colormap jet
title('������� ���: ���������')
xlabel('fast time, us')
ylabel('slow time, s') 
grid on

figure(4)
plot( t_x./1e-6, abs(u(1, :)))
title('�������� ������� ���: ���������')
xlabel('fast time, us')
grid on

figure(5)
plot( t_x./1e-6, real(u(1, :)))
title('�������� ������� ���: real')
xlabel('fast time, us')
grid on

figure(6)
plot( t_x./1e-6, imag(u(1, :)))
title('�������� ������� ���: imag')
xlabel('fast time, us')
grid on
% 
% figure(7)
% plot( t_y, imag(u(:, 4500)))
% title('����������� ������ ���: imag')
% xlabel('fast time, us')
% grid on
  