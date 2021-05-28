%% squint SAR model 24/05/2021
% aleksei.rostov@protonmail.com
clc
clear
close all

fprintf(">> Start SAR Model \n");
%% ���������
gr       = 180 / pi;
c        = 3e8;
%% ��������� ����
Vsar     = 250;   % �������� �������� ����
zsar     = 10000; % ������ ������
% ����� ������� �����������������
x0       = 0000; 
y0       = 8000;
z0       = 0;

fc       = 10e9;        % ������� ������������ �������
Lam      = c/fc;        % ����� �����
La       = 0.8;         % ������ �������� �������� �������� ���
Teta05   = Lam/La*gr;   % ������ ��� �� ������ 0.5
fprintf(">> ������ �� ���� ��� %2.2f ���� \n", Teta05);

dl       = 1; % �, ��������� ����������  �� �������
dr       = 1; % �, ��������� ����������  �� ���������

R0       = sqrt(x0^2 + y0^2 + zsar^2);
fprintf(">> ��������� ��������� ������ ������� ����������������� %5.2f � \n", R0);
TetaQ    = atan(y0/x0)*gr; %
sinTeta0 = (x0^2 + y0^2)/y0^2;
fprintf(">> ���� ������� �� ���� ��� %2.2f ����\n", TetaQ);
%% ����������� ������


Tsyn = Lam*R0 / (dl*2*Vsar*sinTeta0);
fprintf(">> ����� �������������� %1.2f � \n", Tsyn);

dev  = c/(2*dr);  % �������� ���  ��������
tau  = 1e-6;      % ������������ ��� ��������
dt   = 1/dev/1;   % ������ ������ ���
dxI  = dr;        % ��� �� ���������
dyI  = dl;        % ��� �� �������


%% ������������ ������� ������� ����������
% ���������� ����� ���������� � ������ ������������ ������ �������
% ����������������� x0 y0 z0
% Ntarget   - ���������� ��������� ����������
% Map_xyzF - ���������� ��������� � �� ���:     Map_xyzF(1, :) - x axis
%                                               Map_xyzF(2, :) - y axis
%                                               Map_xyzF(3, :) - z axis
%                                               Map_xyzF(4, :) - F ���
Ntarget   = 1;
Map_xyzF = zeros(4, Ntarget);
xi = [20  0 10];
yi = [0 -10 10];
zi = [0   0  0];
Fi = [1   1  1];

Ntarget = length(xi);

%% TODO
% ���������� ������ ����������
% ��� �������� ������� ������������ ������� 
% � ������ ���� ���������
Fprf = 700;
Tp   = 1/Fprf;

if (Fprf > c/(2*R0))
    fprintf(">> �������� �� ��������� \n");
    return
elseif(Fprf < 2*Vsar/La)
    fprintf(">> �������� �� ������� \n");  
    return
end

fprintf(">> ����� ������� ���������� Fprf\n");
fprintf(">> %8.2f <= %d <= %8.2f \n", 2*Vsar/La, Fprf, c/(2*R0));

%% ������������ ��������� �� ������� � ���������
% ������ ������� ��������������
Ls   = R0*(tan((90 - TetaQ+.5*Teta05)/gr) - tan((90 - TetaQ-.5*Teta05)/gr));% length aperture
fprintf(">> ������ ������� �������������� %5.2f � \n", Ls);
% ������������ ��������� ��������������
% Tsyn = Ls/Vsar;
fprintf(">> ������������ ��������� �������������� %2.2f � \n", Tsyn);
% ����������� ����������� ���������
My   =   ceil(Tsyn/Tp);
% ��������� "����������" �������
ty   = (0 : My-1)*Tp;
% ����� ������ � ��������� ������ ����������� �������
tmin = 2*R0/c - 1.5*tau;
tmax = 2*R0/c + 1.5*tau;
% ���������� �������� �� ���������
Mx   = 2*ceil((tmax - tmin)/2/dt);
% ��������� "��������" �������
tx   = tmin + (0:Mx-1)*dt;
fprintf(">> ������� ��� Mx = %d, My = %d \n", Mx, My);

%% ������������ ����������� ������� (������� ���)
% ������� ���
s_raw = zeros(My, Mx);
tn    = zeros(1, My);

for ny = 1 : My
      tn(ny) = ty(ny) - Tsyn / 2;
    for m = 1 : Ntarget 
          R  = sqrt((x0 + xi(m) - Vsar*tn(ny))^2 + (y0 + yi(m))^2 + zsar^2 );
          td = tx - 2*R/c; 
s_raw(ny, :) = s_raw(ny, :) + Fi(m).*exp(1i*pi*dev/tau*(td.^2-td*tau))*exp(-1i*4*pi*R./Lam).*(td>=0 & td<=tau);
    end
    if(mod(ny, 100) == 0)
    fprintf("...%d", round(ny/My*100));
    end
end

fprintf("\n");

figure
imagesc(real(s_raw))
title("SAR Raw signal: real part")
xlabel('Range bins')
ylabel('Azimuth bins')
grid on


%% conv in Range direction
td0      = tx - 2*R0/c;
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
imagesc(tx.*c/2, 1:My, abs(s_range))
title('Range compression')
xlabel('range time bins')
ylabel('azimuth time bins')
grid on

%% range cells correction
% FFT length
NAzFFT = My;
fa     = (1 : NAzFFT);
% ��� �� �������
fsmb  = zeros(NAzFFT, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l), NAzFFT)); 
end

%% conv in Azimuth direction
% �������� ������������ �������
Ka    = 2*Vsar^2/(Lam*R0);
% ���� ������������ �������
FazOp   =  1i*pi*Ka.*ty.*(2*ty(round(My/2+1))-ty);
smb0    =  exp(FazOp);
fsmb0   =  fftshift(fft(smb0, NAzFFT)).'; %     
fsac    =  zeros(NAzFFT, Mx); 
sac     =  zeros(NAzFFT, Mx); 
sF_range = zeros(NAzFFT, Mx); 
for l = 1 : Mx
    fsac(:,l) =  fsmb(:, l).*conj(fsmb0);  
     sac(:,l) =  fftshift(ifft(fsac(:,l))); 
end

figure
imagesc(tx.*c/2, tn.*Vsar,abs(sac))%  
title('Radar Image')
xlabel('range: m')
ylabel('azimuth: m')
grid on


fprintf(">> End SAR Model \n");


