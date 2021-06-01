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
x0       = 4000; 
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

Tsyn = Lam*R0 / (dl*2*Vsar*sinTeta0) + 0.0;
fprintf(">> ����� �������������� %1.2f � \n", Tsyn);

dev  = c/(2*dr);  % �������� ���  ��������
tau  = 1e-6;      % ������������ ��� ��������
dt   = 1/dev/1;   % ������ ������ ���
dxI  = dr;        % ��� �� ���������
dyI  = dl;        % ��� �� �������


%% ������������ ������� ������� ����������
% ���������� ����� ���������� � ������ ������������ ������ �������
% ����������������� {x0 y0 z0}
% Ntarget   - ���������� ��������� ����������
% Map_xyzF - ���������� ��������� � �� ���:     Map_xyzF{1} - x axis
%                                               Map_xyzF{2} - y axis
%                                               Map_xyzF{3} - z axis
%                                               Map_xyzF{4} - F ���
xi = [20  0];
yi = [20  0];
zi = [ 0  0];
Fi = [ 1  1];
% 

% xi = [0];
% yi = [0];
% zi = [0];
% Fi = [1];

Map_xyzF = cell(4, 1);
Ntarget  = length(xi);

if(Ntarget ~= (length(yi) + length(zi) + length(Fi))/3 )
    fprintf(">> �� ��������� ������ ���� \n");
    return
end

Map_xyzF{1} = xi; 
Map_xyzF{2} = yi;
Map_xyzF{3} = zi;
Map_xyzF{4} = Fi;



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
Tsyn = Ls/Vsar;
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
Rref = zeros(1, My);
for ny = 1 : My
      tn(ny) = ty(ny) - Tsyn / 2;
     Rref(ny)  = sqrt((x0 - Vsar*tn(ny))^2 + (y0)^2 + zsar^2 );

    for m = 1 : Ntarget 
          R  = sqrt((x0 + Map_xyzF{1}(m) - Vsar*tn(ny))^2 + (y0 + Map_xyzF{2}(m))^2 + zsar^2 );
          td = tx - 2*R/c; 
s_raw(ny, :) = s_raw(ny, :) + Map_xyzF{4}(m).*exp(1i*pi*dev/tau*(td.^2-td*tau))*exp(-1i*4*pi*R./Lam).*(td>=0 & td<=tau);
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
     fs_raw(k , :)   = fft(s_raw(k, :).*hann(Mx).');
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

%  fa     = (1 : NAzFFT);

fa = -1/Tp/2:1/Tsyn:1/Tp/2;
dD = R0*(sqrt(fa.^2*Lam^2/(4*Vsar^2)+1)-1);

rangD    = round(dD/dxI);
rangDmax = max(rangD);

range_dR = round((Rref - min(Rref))/dxI);

% ������� ��� ��������� 
% N_1 = R3*Lam^2/(4*Vsar^2/dxI);
% dNr = ceil(N_1.*(fD1 - fa./(Tp*NAzFFT)).^2 - N_1*fD1^2);
Tsmb  = zeros(NAzFFT, Mx);
for k = 1 : My
        Tsmb(k,:) = circshift(s_range(k,:),  -range_dR(k));
end


% ��� �� �������
fsmb  = zeros(NAzFFT, Mx);
for l=1:Mx
    fsmb(:,l)=fftshift(fft(Tsmb(:,l), NAzFFT)); 
end
fsmbF = fsmb;


figure
imagesc(abs(fsmb))
title('before RCMC: FFT azimuth')
grid on

 
for k = 1 : My
        fsmbF(k,:) = circshift(fsmb(k,:), -range_dR(k));
end

%  
% for k = 1 : My
%     for m = 1 : Mx - rangDmax 
%         fsmbF(k,m) = fsmb(k, m + rangD(k));
%     end
% end

figure
imagesc(abs(fsmbF))
title('after RCMC: FFT azimuth')
grid on


smbF = zeros(NAzFFT, Mx);
for l = 1 : Mx
    smbF(:,l)=ifft(fftshift(fsmbF(:,l)));
end


figure
imagesc(tx.*c/2, 1:My, abs(smbF))
title('Range cells migration correction')
xlabel('range time bins')
ylabel('azimuth time bins')
grid on

%% conv in Azimuth direction
% �������� ������������ �������
Ka    = 2*Vsar^2/(Lam*R0);
% ���� ������������ �������
FazOp   =  1i*pi*Ka.*ty.*(2*ty(round(My/2+1))-ty);

smb00 = exp(-1i*4*pi*Rref./Lam);

winvec = hilbert(hann(My));

smb0    =  exp(FazOp);
fsmb0   =  fftshift(fft(smb00, NAzFFT)).'; %     
fsac    =  zeros(NAzFFT, Mx); 
sac     =  zeros(NAzFFT, Mx); 
sF_range = zeros(NAzFFT, Mx); 
for l = 1 : Mx
    fsac(:,l) =  fsmb(:, l).*conj(fsmb0);  
     sac(:,l) =  fftshift(ifft(fsac(:,l))).*winvec; 
end

figure
imagesc(tx.*c/2, tn.*Vsar,abs(sac))%  
title('Radar Image')
xlabel('range: m')
ylabel('azimuth: m')
grid on

figure
mesh(abs(sac))


fprintf(">> End SAR Model \n");


