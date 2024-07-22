function [L,Gt] = logLikelihood(thetaS, thetaD, outputs, model, specification)
    % Tamaño de los datos
   
    [i,j]=models_index(model, specification)
    outputs
    rotated_returns=outputs.rotated_returns
    T = size(rotated_returns, 1);
    d = size(rotated_returns, 2);
    P = outputs.P
    Lambda=outputs.Lambda
    Dt= outputs.Dt

    % Prealocar Gt, Qt y Ct para eficiencia
    Gt = zeros(d, d, T+1); % T+1 para inicialización en t=0
    Gt(:, :, 1) = eye(d); % Valor inicial de Gt para t=0 %abajo sigue la iteraci�n hasta completar Gt
    Id = eye(d); % Identidad de tamaño d
    
    if i==3

    % Calcular Qt_star, Qt y Ct para la log-verosimilitud de RDCC
    Qt_star = calc_all_Gts(thetaS, thetaD, rotated_returns, Id, model, specification);
    Qt = zeros(d, d, T+1);
    Ct = zeros(d, d, T+1);
    
    Qt(:, :, 1) = P * sqrt(Lambda) * P' * Qt_star(:, :, 1) * P * sqrt(Lambda) * P';
    Ct(:, :, 1) = sqrt(inv(Qt(:, :, 1) .* Id)) * Qt(:, :, 1) * sqrt(inv(Qt(:, :, 1) .* Id));
    end
    
    for t = 2:T+1
        Qt(:, :, t) = P * sqrt(Lambda) * P' * Qt_star(:, :, t) * P * sqrt(Lambda) * P';
        Ct(:, :, t) = sqrt(inv(Qt(:, :, t) .* Id)) * Qt(:, :, t) * sqrt(inv(Qt(:, :, t) .* Id));
    end

    % Inicializar la log-verosimilitud
    ll = 0;
    L = ll_type(model, d, Gt(:, :, 1), rotated_returns, 2,Dt,Ct);
    
    for t = 2:T+1 % Ajustando el bucle para que t comience desde 2
        Gt(:, :, t) = calcGt(thetaS, thetaD, rotated_returns(t-1, :), Gt(:, :, t-1), model, specification); % t-1 para los retornos para coincidir con la teoría
        
        % Agregar un pequeño término regularizador a Gt para evitar singularidades
        reg_term = 1e-6 * eye(d);
        Gt_reg = Gt(:, :, t) + reg_term;

        ll = ll_type(model, d, thetaS, Gt_reg, rotated_returns, t)
        disp('estoy en logLikelihood, volvele a poner el punto y coma');

        L = L + ll;
    end
end
