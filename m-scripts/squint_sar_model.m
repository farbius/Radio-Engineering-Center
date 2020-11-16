%%
%
clc
clear
close all
LOAD = 1;
fprintf(">> Start SAR Model \n");
%% ���������
gr   = 180 / pi;
c    = 3e8;
%% ��������� ����
Vsar = 250;   % �������� �������� ����
zsar = 10000; % ������ ������
% ����� ������� �����������������
x0   = 10000; 
y0   = 0;
z0   = 0;

fc    = 10e9;        % ������� ������������ �������
Lam   = c/fc;        % ����� �����
La    = 1.0;         % ������ �������� �������� �������� ���
Teta05= Lam/La*gr;   % ������ ��� �� ������ 0.5
fprintf(">> ������ �� ���� ��� %2.2f ���� \n", Teta05);

TetaQ = 5; %
fprintf(">> ���� ������� �� ���� ��� %2.2f ����\n", TetaQ);

%% ����������� ������
dev  = 150e6;  % �������� ���  ��������
tau  = 6e-6;   % ������������ ��� ��������
dt   = 1/dev/2;% ������ ������ ���
dxI  = c/dev;  % ��� �� ���������
dyI  = La/2;   % ��� �� �������



%% ������������ ������� ������� ����������
targets = rgb2gray(imread('targets', 'png'));
[Nx, Ny] = size(targets);
Ntarget = Nx*Ny;
num=1; xn=zeros(Ntarget,1); 
yn=xn; Fn=xn; 
% ������������ �������� ��������� ������� �����������������
% ���������� z ����������� 
for iNy=1:Ny 
    for iNx=1:Nx 
        xn(num)=(iNx-Nx/2);   
        yn(num)=(Ny/2-iNy+1);  
        Fn(num)=double(targets(iNy,iNx))/1;     
        num=num+1; 
    end 
end 
% ����������� ����������
kx = 3; ky = 3;
xn=xn*kx*dxI; yn=yn*ky*dyI;


%% ��������� ��������� ������ ������� ����������������� � ��������� ����
R_0  = sqrt(x0^2 + zsar^2);
% ��������� ���� �� ��������� ��������������
u1   = -R_0*tan((TetaQ+.5*Teta05)/gr); % ������ ���������
u2   = -R_0*tan((TetaQ)/gr);           % �������� ���������
u3   = -R_0*tan((TetaQ-.5*Teta05)/gr); % ��������� ���������
fprintf(">> ��������� ���� �� ��� u\n");
fprintf(">> u1 = %5.2f �, u2 = %5.2f �, u3 = %5.2f �,\n", u1, u2, u3);
% ��������� ��������� �� ��������� ��������������
R3   = sqrt(u3.^2   + R_0^2);
R2   = sqrt(u2.^2   + R_0^2);
R1   = sqrt(u1.^2   + R_0^2);
fprintf(">> ��������� ��������� R1 = %8.2f m, R2 = %8.2f m, R3 = %8.2f m\n", R1, R2, R3);
% ������ �� ��������� ��������������
fD1  = 2*Vsar/Lam*sin((TetaQ-.5*Teta05)/gr);
fD3  = 2*Vsar/Lam*sin((TetaQ+.5*Teta05)/gr);
fD   = fD3 - fD1;
fD2  = fD1 + fD/2;
fprintf(">> ������ ������ fD1 = %8.2f �� \n", fD1);
fprintf(">> ������� ������ fD3 = %8.2f �� \n", fD3);
fprintf(">> ������� ������ fD2 = %8.2f �� \n", fD2);
fprintf(">> ������ ������� �� ������� %8.2f �� \n", fD);

%% TODO
% ���������� ������ ����������
% ��� �������� ������� ������������ ������� 
% � ������ ���� ���������
Fprf = 700;
Tp   = 1/Fprf;

if (Fprf > c/(2*abs(R1)))
    fprintf(">> �������� �� ��������� \n");
    return
elseif(Fprf < 2*Vsar/La)
    fprintf(">> �������� �� ������� \n");  
    return
end

fprintf(">> ����� ������� ���������� Fprf\n");
fprintf(">> %8.2f <= %d <= %8.2f \n", 2*Vsar/La, Fprf, c/(2*abs(R1)));

%% ������������ ��������� �� ������� � ���������
% ������ ������� ��������������
Ls   = R_0*(tan((TetaQ+.5*Teta05)/gr) - tan((TetaQ-.5*Teta05)/gr));% length aperture
fprintf(">> ������ ������� �������������� %5.2f � \n", Ls);
% ������������ ��������� ��������������
Tsyn = Ls/Vsar;
fprintf(">> ������������ ��������� �������������� %2.2f � \n", Tsyn);
% ����������� ����������� ���������
My   =   ceil(Tsyn/Tp);
% ��������� "����������" �������
ty   = (0 : My-1)*Tp;
% ����� ������ � ��������� ������ ����������� �������
tmin = 2*R_0/c - 0.5*tau;
tmax = 2*R_0/c + 1.5*tau;
% ���������� �������� �� ���������
Mx   = 2*ceil((tmax - tmin)/2/dt);
% ��������� "��������" �������
tx   = tmin + (0:Mx-1)*dt;
fprintf(">> ������� ��� Mx= %d, My= %d \n", Mx, My);
%%
% ��������� ��������� ���������
u = linspace(u1, u3, My);
R_center   = sqrt(u.^2 + R_0^2);
Teta = linspace(TetaQ+.5*Teta05, TetaQ-.5*Teta05, My);
f_dop= 2*Vsar/Lam*sin(Teta./gr);

figure
plot(R_center, f_dop, '.-b', R3, fD1, 'or', R2, fD2, 'or', R1, fD3, 'or')
txt = ['Teta_H: ' num2str(TetaQ) ' ����'];
text(R2+1.0,fD2,txt)
title('��������� ������� �� ��������� ��������������')
xlabel('Range, m')
ylabel('������, ��')
grid on

%% ������������ ����������� ������� (������� ���)
% ���������� u ���� �� ��������� ��������������
u = linspace(u1, u3, My);
% ������� ���
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
% ��� �� �������
fsmb  = zeros(NAzFFT, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(s_range(:,l), NAzFFT)); 
end
% ������� ��� ��������� 
N_1 = R3*Lam^2/(4*Vsar^2/dxI);
dNr = ceil(N_1.*(fD1 - fa./(Tp*NAzFFT)).^2 - N_1*fD1^2);
% broadside view
% dNr = R_0*sqrt(1 - (k.*Lam*Fprf/(2*Vsar*NAzFFT)).^2) - R_0;
% ��������� �������� ���������
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
% �������� ������������ �������
Ka    = 2*Vsar^2/(Lam*R_0);
% ���� ������������ �������
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

