%% ������������ ������� ������������ �����������
% ������ �������, 25.05.20
% farbius@protonmail.com

% �������� �� ����������
x_size = (floor(x_size*.5))*2 + 1;
y_size = (floor(y_size*.5))*2 + 1;
% ����� ��������� ����������
Num_trg   = x_size * y_size;

% ����������� ���������� ������� �������� � ������ ��
% ����� ������� -- ��� ����������� ����� (��) ������������ �����������
x_p   = linspace(-floor(x_size/2), floor(x_size/2), x_size)';
y_p   = linspace(-floor(y_size/2), floor(y_size/2), y_size);
x_p   = repmat(x_p, [1 x_size]);
y_p   = repmat(y_p, [y_size 1]);

% ����������� ���������� ��������� �� ������������ �������� ������ �������
% (���) ��� (�������� � ������� ������ ������� ��� ����� ���������� (0,0,z0))
Rn    = sqrt(R_0^2 - z0^2); % �������� ��������� ��������� � ��������
% ���������� �� � ������ ������� ��������� ������������ ���
x_cp  = sqrt(R_0^2 - z0^2)*cos(PsiAz*gr);
y_cp  = sqrt(R_0^2 - z0^2)*sin(PsiAz*gr);
% ��������� ������������ ����������� � ������ ��������� ��� � ������ ��
x_p = x_p.*dx + x_cp;
y_p = y_p.*dy + y_cp;

y_p = y_p(:);
x_p = x_p(:);

% �������� ������ �������� �������� (���) ��������� ����������
% � ������ ������� ������ ��������� ��� ���� �������� ����������
% ��� ������� ������������� ����� �������� ������ ���� My - �� ����������
% ������������
Sigma = zeros(x_size,y_size);
% ������ ����������� ������� �������
Sigma(round(x_size/2) + 1, round(y_size/2) + 1) = 1;
Sigma = Sigma(:);

