function [Qt, Qt_star , Ct] = calcQt(model,specification, outputs, thetaD)

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
    Id=eye(d);

    % Prealocar Qt y Ct para eficiencia
        
    % Inicializar Qt y Qt_star si i==3
    if i == 4 % 4 is the model index for RDCC
        % Calcular Qt y Ct para t = 1 (sería t=0 en teoría)
        initial_Qt_star=Id
        Qt_star = calc_all_Gts(model,specification,outputs, thetaD, initial_Qt_star); uso la función de G porque es la misma 
        Qt = zeros(d, d, T);
        initial_Qt = P * sqrt(Lambda) * P' * initial_Qt_star * P * sqrt(Lambda) * P';
        
        Ct = zeros(d, d, T);
        initial_Ct = sqrt(initial_Qt .* Id) * initial_Qt * sqrt(inv(initial_Qt .* Id));

        % Calcular Qt y Ct para t >= 1
        for     t=1:T
        Qt(:, :, t) = P * sqrt(Lambda) * P' * Qt_star(:, :, t) * P * sqrt(Lambda) * P';
        Ct(:, :, t) = sqrt(inv(Qt(:, :, t) .* Id)) * Qt(:, :, t) * sqrt(inv(Qt(:, :, t) .* Id));
        end
    else
        Qt = [];
        Qt_star = [];
        Ct=[]
    end