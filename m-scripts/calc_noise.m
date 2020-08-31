%% формирование матрицы внутренних шумов
% Ростов Алексей, 25.05.20
% farbius@protonmail.com

n = zeros(My, Mx);
for ny = 1 : My
    for m = 1 : Num_trg
        n(ny,:) = n(ny,:) + (randn(1,Mx) + 1i*randn(1,Mx)).*std_noise;
    end
end