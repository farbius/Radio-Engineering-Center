%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                    % 
% ������ ����������������� ��� JSTARS � �������� �����������   %      
% ����������� �����                                            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clc; clear all;
gr=pi/180; c=3e8;
 
CALC_SIGNAL=0;  % 1 - ������������ ������; 
                % 0 - ������������ ����������� � ����� work_s.mat
                % ------ ����. �������  <s, n, F, tn>  -------
Signal=1;       % 1 - �������� ������ ����; 0 - ���
Noise=0;        % 1 - ���������� ��� ����; 0 - ���
Jamm=1;         % 1 - ������ ����; 0 - ���
 
%��������� ���
Lam=0.03; %�, ����� �����
Lam2_dB=20*log10(Lam);
fo=c/Lam; %��, ������� �������
Tet05=1.6;  %��, ������ �������� ��
% Eps=10; %��, ���� ������� �� � ����. ���������
dev=150e6; %��, �������� ���
tau=4.0e-6; %c, ������������ �������� ���    
Tp=2e-3; %c, ������ ���������� ���������      
Vn=200; %�/�, �������� �������� ��������
D0=60e3; %�, �������������� ��������� ������ ������� ����������������� �� ��������
La=Lam/(Tet05*gr); %�, ���������� ������ �������� �������
G=10000; %�� �������
G_dB=10*log10(G);
Pizl=10000; %��, ���������� �������� ���������
Pizl_dB=10*log10(Pizl);
k4pi_dB=10*log10(4*pi);
kn_dB=6;%�-� ���� ��������� ����������
SigZ_dB=10; %,dB ��.�., ��� ���� (����. ��� ������������ ��������)
Sig0_dB=-20; %,dB ��.�., ��� ������� ������ ����������� � �������� ����������
 
%��������� ������� �����
nJ=103; %����� �������� ��������� �����������, � ������� ����������� ������� �����   6*16+7=103 
% Tet05J=60; %��, ������ �������� �� �������� � ���������� ������
Tet05J=1.6; %��, ������ �������� �� �������� � ���������� ������ ���-30
KuJ_dB=100; %��, �-� �������� ������ ������������
PmaxJ=20; %��, ������������ �������� �������� ���
Pprd_spn=640; % ��, �������� ��� ���-30
Pprd_spn_dB=26;%10*log10(Pprd_spn);
Fsh=30e6; % ���, ������ ������� ������� ������ ���-30
Fsh_dB=10*log10(Fsh); % �����, ������ ������� ������� ������ ���-30
 
TypeJD=0 %��� ������ ������ ���������:
         %              0 - ��� ������
         %              1 - ���������� ����� �� ������� �� �������� delFJ
         %              2 - ??? ���������� ����� �� ������� �� �������� delFJ
         %                  �� ��������� ���. ����� �� ������������ � ������������
         %              3 - 0/pi ���������� ���� � �������� T0piJ
         %              4 - 0/pi ���������� ���� �� ������
         %                  �-������������������ � ���������� L_mp ��������� �
         %                  ������� ������������ ��������
         %              5 - ������� ������ �� ���-30
 
TypeJA= 2 %��� ������ ������ �������:
         %              0 - ��� ������
         %              1 - ���������� ����� � ������������ ������ ��
         %              �������� delFJ_Az
         %              2 - 0/pi � ������������ ������ � �������� T0piJA
         %             
delFJ=-0e6; %���������� ����� �� ������� (TypeJD=1,2) 15 ���
T0piJ=0.2e-6; %������������ ����. �������� 0/pi ���������� (TypeJD=3)
L_MP=10; %��������� M-������������������   (TypeJD=4)
delFJA=-15.0; %��, ��������� ����� � ������������ ������ (TypeJA=1)
T0piJA=0.2; %������������ ����. �������� 0/pi ���������� (TypeJA=2)
 
GJ=36000/Tet05J^2; %�� ������� 
GJ_dB=10*log10(GJ);
PmaxJ_dB=10*log10(PmaxJ);
 
%��������� ������������ �����������
dxI=1/dev*150e6/2; %��� ����� ��� ��������� - �������� ����������� ����������� �� D
dyI=La/2; %��� ����� ��� ������� - �������� ��������
 
% ������������� ��������� ����������� 
target_name='pole11'; %�������� ����� ��������� ����������������� �������  (GIF Grayscale Image File) 
target=imread(target_name,'gif'); % ������ ������� ����������� 
[Nx Ny]=size(target); %����������� �� ����
Ntarget=Nx*Ny;        %���������� ����������� ��������� 
kx=12; %�-� ���������� ��������� ����������� ����� ����� ��������� �� ��������� � ���� dxI 
ky=12; %�-� ���������� ��������� ����������� ����� ����� ������� �� ��������� � ���� dyI 
 
Ls=1600;%D0*Tet05*gr; %�, ����� ������� �������������� (����� ����� �������)
Tsint=Ls/Vn; %c, ����� ��������������;
 
Mxi=kx*(Nx-1)+1; % ���-�� ��������� ���������� �� D, ���������� ������������
Myi=ky*(Ny-1)+1; % ���-�� ��������� ���������� �� ��, ���������� ������������
My=floor(Tsint/Tp); % ���-�� ��������� ���������� �� �� (����� ������������)
 
num=1; xn=zeros(Ntarget,1); yn=xn; Fn=xn; % ������������� ���������
for iNy=1:Ny 
    for iNx=1:Nx 
        xn(num)=(iNx-Nx/2);   % ������ ������ �� ���������� x �������� (������)
        yn(num)=(Ny/2-iNy+1); % ������ ������ �� ���������� y �������� (������� �� ������ �� ������)
        fn(num)=double(target(iNy,iNx))/251; % ��������� ��������� �������� �� �������
        num=num+1; 
       
    end 
end 
xn=xn*kx*dxI; yn=yn*ky*dyI; % ��������� � ������ ����������
 
Ts=(2*(D0+min(xn)))/c - 3.5e-6; % ��������� ����� ������� �� ����� ���������
Tf=(2*(D0+max(xn)))/c + 7.5e-6; %tau ��������  ����� ������� �� ����� ���������
dt=1/dev/1; %��� ������� �� ������� - �������� ����. ����������� ������� �� ���������
Mx=2*ceil((.5*(Tf-Ts))/dt); %���������� ������� �� ����� ���������
t= Ts+(0:Mx-1)*dt; % ������ ������� ������� �� ����� ���������
ty=[0:My-1]*Tp; %������ ������� �� �������������
Fd=dev; % ���, ������� �������������
Fd_dB=10*log10(Fd); % �����
 
Kr=dev/tau; %�������� ���
Ka=2*Vn^2/(Lam*D0); % �������� ��� ��� ������������ ���������
 
%������ �������� ����
dFn_MHz_dB=10*log10(1/(tau*My)/1e6); %������� ������ ���������, ����� 
Pn_dB=-144+kn_dB+dFn_MHz_dB; %�������� �����. ���� ���������
 
%������ �������� ������������ ������� (��� �������� � ��� 1 ��.�)
D0_dB=10*log10(D0); %���������, ��
P1m2_dB=Pizl_dB+2*G_dB+Lam2_dB-3*k4pi_dB-4*D0_dB;
q_dB=P1m2_dB-Pn_dB %��������� ������-���  ��� �������� � ��� 1 ��.�
q=10^(q_dB/20);%��������� ������-��� �� ��������� ��� �������� � ��� 1 ��.�
 
%������� ����������� ����������� � ��. ���
SigXY_dB=fn*(SigZ_dB-Sig0_dB)+Sig0_dB; 
SigXY=10.^(SigXY_dB/10);
 
% ������������ ����������� ������ ����������� �������, ���� � ������
s=zeros(My,Mx);   n=zeros(My,Mx); jm=zeros(My,Mx); u=zeros(My,Mx); 
sqrt2=sqrt(2);
Xground=zeros(My,Mx); 
 
if CALC_SIGNAL==0
    load('work_sss.mat');
else
    for ny=1:My %������ �� ������ ������������
        tn(ny)=ty(ny)-Tsint/2;
        F(ny)=sinc(La*(atan(Vn*tn(ny)/D0))/Lam); %sin(x)/x �� �������� �������
        for m=1:Ntarget % ������ �� ������ ������������ �������
            tnm=tn(ny)+yn(m)/Vn; %����� ��������� m-�� ������� � n-� ������������
            R=sqrt((D0+xn(m))^2+Vn^2*(tnm^2)); %���������� �� �������. ������ ������ �������� �� ������!!!!!
            td=t-2*R/c; %������ �������, ��������� � �����. �������. td>0 � ������ ������� ������ ���. �������
            %����� ��������, ���������� �� �������� � ����� ������������
            s(ny,:)=s(ny,:)+q*SigXY(m)*F(ny).^2*exp(-j*(4*pi*R/Lam-pi*Kr*(td.^2-td*tau))).*(td>=0 & td<=tau);
            %���
            n(ny,:)=n(ny,:)+(randn(1,Mx)+j*randn(1,Mx))/sqrt2;
            Xground(1,:) = fn(m);
         end; 
         if mod(ny,10)==0 
             display(ny/My*100); 
         end;  % m=1:Ntarget
    end;  % ny=1:My
end;  % CALC_SIGNAL
 
if Signal==1 
    u=s; 
end;
 
if Noise==1 
    u=u+n; 
end;
 
 if Jamm==1
 %������ �������� ������
F_dB=20*log10(F); %�� ������� ��� �� �������������
PprmJ_dB=Pizl_dB+G_dB+F_dB-2*D0_dB+GJ_dB+Lam2_dB-2*k4pi_dB; %�������� ������� �� ������ �������� �������
PprdJ_dB=PprmJ_dB+KuJ_dB;
PprdJ0_dB=PprdJ_dB;
 
NM=2^L_MP-1;
tauM=tau/NM; %������������ ������������ ������� �-������������������
MP=M_POSL(L_MP,0); %����� ��������� ���� ������ � ������� ������������ ����. ��������
 
for ny=1:My %������ �� ������ ������������
    if PprdJ_dB(ny)>PmaxJ_dB
        PprdJ_dB(ny)=PmaxJ_dB; %����������� ������������ �������� ����������� ������
    end;
end;
PJ_dB=PprdJ_dB-2*D0_dB+G_dB+F_dB+Lam2_dB-2*k4pi_dB; %�������� ������ � �������� ������
 
qJ_dB=PJ_dB-Pn_dB; %��, ��������� ������-��� � �������� ������
qJ_dB_max=max(qJ_dB)
qJ=10.^(qJ_dB/20); %
PJ_dB_s=Pprd_spn_dB-2*D0_dB+GJ_dB+G_dB+Lam2_dB-2*k4pi_dB+(Fd_dB-Fsh_dB); %�������� ������ � �������� ������
% ��������� ������-���
qJ_dB_spn=PJ_dB_s-Pn_dB;
q_SP_dB=q_dB-qJ_dB_spn;
qJ_spn=10.^(qJ_dB_spn/20);
 
    for ny=1:My %������ �� ������ ������������
        tnm=tn(ny)+yn(nJ)/Vn; %����� ��������� �������� 
        RJ=sqrt((D0+xn(nJ))^2+Vn^2*(tnm^2)); %���������� 
        td=t-2*RJ/c; %������ �������, ��������� � �����. �������. td>0 � ������ ������� ������ ���. �������
        % ������ �������� ������� ������  �� ������ ��� �������
        
        
        %����������������� ������ ��� ���������
        jm(ny,:)=qJ(ny)*exp(-j*(4*pi*RJ/Lam-pi*Kr*(td.^2-td*tau))).*(td>=0 & td<=tau);
        
        if TypeJD==1 %���������� ��������� ��������
            jm(ny,:)=jm(ny,:).*exp(j*2*pi*delFJ*td).*(td>=0 & td<=tau);
        end; %TypeJD=1
 
        if TypeJD==2 %���������� ��������� �������� �� ����. ���. �����
            psi=rand*2*pi;
            jm(ny,:)=jm(ny,:).*exp(j*(2*pi*delFJ*td+psi)).*(td>=0 & td<=tau);
        end; %TypeJD=2
        
        if TypeJD==3 %0/pi ���������� ����
            Ni=floor(Mx*dt/T0piJ); %���-�� ����� �������� ���������� ���� �� ������������ �������
            Ki=floor(Mx/Ni); %���������� ����� �������� ������� � 1 ������� ���������� ����
            FAZ_0pi=zeros(1,Mx);     FAZ_0pi(1)=1;
            for ki=2:Mx
                if mod(ki,Ki)==0
                    FAZ_0pi(ki)=-FAZ_0pi(ki-1);
                else
                    FAZ_0pi(ki)=FAZ_0pi(ki-1);
                end;
            end;
            jm(ny,:)=jm(ny,:).*FAZ_0pi.*(td>=0 & td<=tau);
        end; %TypeJD=3
        
        if TypeJD==4 %0/pi ���������� ���� �� ������ �-������������������
            Ni=NM; %���-�� ����� �������� ���������� ���� �� ������������ �������
            Ki=floor(Mx/Ni); %���������� ����� �������� ������� � 1 ������� ���������� ����
            FAZ_0pi=zeros(1,Mx);   iM=1;  FAZ_0pi(1)=MP(iM);
            
            for ki=2:Mx
                if mod(ki,Ki)==0
                    if iM<NM
                        iM=iM+1;
                    end;
                    FAZ_0pi(ki)=MP(iM);
 
                else
                    FAZ_0pi(ki)=FAZ_0pi(ki-1);
                end;
            end;
            jm(ny,:)=jm(ny,:).*FAZ_0pi.*(td>=0 & td<=tau);
        end; %TypeJD=4
        
         if TypeJD==5 % ������� ������ �� ���-30
            
           jm(ny,:)=qJ_spn.*(randn(1,Mx)+j*randn(1,Mx))/sqrt2;
            
        end; %TypeJD=5
 
        if TypeJA==1 %������������ �������� ����� 
            jm(ny,:)=jm(ny,:).*exp(j*2*pi*delFJA*tnm).*(td>=0 & td<=tau);
        end; %TypeJA==1
       
        if TypeJA==2 %0/pi ���������� ���� 
            Ni=floor(Tsint/T0piJA); %���-�� ����� �������� ���������� ���� �� ������������ c�������������
            Ki=floor(My/Ni); %���������� ����� ������������ � 1 ������� ���������� ����
            FAZ_0piA=zeros(1,My);     FAZ_0piA(1)=1;
            for ki=2:My
                if mod(ki,Ki)==0
                    FAZ_0piA(ki)=-FAZ_0piA(ki-1);
                else
                    FAZ_0piA(ki)=FAZ_0piA(ki-1);
                end;
            end;
            jm(ny,:)=jm(ny,:)*FAZ_0piA(ny).*(td>=0 & td<=tau);
        end; %TypeJA=2
        
         if TypeJA==3 %0/pi ���������� ���� �� ������ �-������������������
            Ni=floor(Tsint/T0piJA); %���-�� ����� �������� ���������� ���� �� ������������ c�������������
            Ki=floor(My/Ni); %���������� ����� ������������ � 1 ������� ���������� ����
            FAZ_0piA=zeros(1,My);   iM=1;  FAZ_0piA(1)=MP(iM);
            for ki=2:My
                if mod(ki,Ki)==0
                    if iM<NM
                        iM=iM+1;
                    end;
                    FAZ_0piA(ki)=MP(iM);
 
                else
                    FAZ_0piA(ki)=FAZ_0piA(ki-1);
                end;
            end;
            jm(ny,:)=jm(ny,:)*FAZ_0piA(ny).*(td>=0 & td<=tau);
        end; %TypeJA=3
 
 
    end;  % ny=1:My
    u=u+jm;
    
end;  % Jamm==1
 
figure(10) 
imagesc(t*150e6/1000,ty*Vn/1000, real(u));
colormap(gray)
ylabel('����� �������, ��')
xlabel('����� ���������, ��')
title('���������� ������')
 
% C����� �� ��������� 
td0=t-2*D0/c; 
pha20=pi*Kr*((td0.^2)-td0*tau); 
s0=exp(j*pha20).*(td0>=0 & td0<=tau); %������� ������ ��� ������������ ����������
fs0=fft(s0);
for k=1:My; 
    fu(k,:)=fft(u(k,:));
    fsm(k,:)=fu(k,:).*conj(fs0);% ���������� � ��������� �������
    smb(k,:)=fftshift(ifft(fsm(k,:)));
end; 
 
figure(20) 
imagesc(t*150e6/1000,ty*Vn/1000,-abs(smb))
colormap(gray)
ylabel('����� �������, ��')
xlabel('����� ���������, ��')
title('������, ������ �� ���������')
 
% �����������
fa=[-1/Tp/2+1/Tsint:1/Tsint:1/Tp/2]; %����� ������������ ������ (�� ����� ����� ������������ �������� ��� ���������� ����������������� �������)
dD=D0*(sqrt(fa.^2*Lam^2/(4*Vn^2)+1)-1);
 
rangD=round(dD/dxI);%������ � ���-�� ��������� ���������� �� D
rangDmax=max(rangD);
 
for l=1:Mx; 
    fsmb(:,l)=fftshift(fft(smb(:,l))); % ���������������� ������������ ������ �������, ������� �� ���������
end; 
 fsmbF=fsmb;
 
 
for k=1:My; 
    for m=1:Mx-rangDmax 
        fsmbF(k,m)=fsmb(k,m+rangD(k)); %C���� ����� �� rangD(k) �����
    end
end; 
 
 
for l=1:Mx; 
    smbF(:,l)=ifft(fftshift(fsmbF(:,l))); % ��������������� ���������������� ������, ������ �� ���������
end; 
 
figure(31) 
imagesc(t*150e6/1000,ty*Vn/1000,-abs(smbF)) 
colormap(gray)
ylabel('����� �������, ��')
xlabel('����� ���������, ��')
title('������� ��� ����� �� �� ���������')
 
 
% ������ �� �������
smb0=exp(j*pi*Ka*ty.*(2*ty(round(My/2+1))-ty)); %������� ����������� ������ ��� ������������ ���������
fsmb0=fftshift(fft(smb0)); % C����� �������� ������������ ������� 
 
for l=1:Mx; 
   fsac(:,l)=fsmb(:,l).*(fsmb0)'; %  ������������ ����������
   sac(:,l)=fftshift(ifft(fftshift(fsac(:,l)))); % �����������
end; 
 
smax = max(max(abs(sac)));
sacc = zeros(My,Mx);
 
for l=1:Mx; 
    sacc(:,l)=abs(sac(:,l))./max(max(abs(sac))); % ���������������� ������������ ������ �������, ������� �� ���������
end; 
 
x = 100;%round(0/10);
y = 100;%round(0/10);
figure(41) 
imagesc(t*150e6/1000,ty*Vn/1000,-abs(sac)) 
colormap(gray)
ylabel('����� �������, ��')
xlabel('����� ���������, ��')
title('����������� ���')