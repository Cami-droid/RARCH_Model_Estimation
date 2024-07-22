function [Qt, Qt_star] = calcQt(model, specification, outputs, thetaD_opt)
     % Inicialización de variables
    [i, j] = models_index(model, specification);
    rotated_returns = outputs.rotated_returns;
    T = outputs.T;
    d = outputs.d;
    thetaS=outputs.H_bar;
    P = outputs.P;
    Gt =outputs.Gt;
    Lambda = outputs.Lambda;
    Dt = outputs.Dt;

    % Prealocar Qt y Ct para eficiencia
        
    % Inicializar Qt y Qt_star si i==3
    if i == 3

        % Calcular Qt y Ct para t = 1 (sería t=0 en teoría)
        Qt_star = calc_all_Gts(model, specificacion, outputs, thetaD, Id);
        Qt = zeros(d, d, T+1);
        Qt(:, :, 1) = P * sqrt(Lambda) * P' * Qt_star(:, :, 1) * P * sqrt(Lambda) * P';
        
        Ct = zeros(d, d, T+1);
        Ct(:, :, 1) = sqrt(inv(Qt(:, :, 1) .* Id)) * Qt(:, :, 1) * sqrt(inv(Qt(:, :, 1) .* Id));

        % Calcular Qt y Ct para t >= 2
        t_count=2:T+1
        Qt(:, :, t_count) = P * sqrt(Lambda) * P' * Qt_star(:, :, t_count) * P * sqrt(Lambda) * P';
        Ct(:, :, t_count) = sqrt(inv(Qt(:, :, t_count) .* Id)) * Qt(:, :, t_count) * sqrt(inv(Qt(:, :, t_count) .* Id));
    else
        Qt = [];
        Qt_star = [];
    end