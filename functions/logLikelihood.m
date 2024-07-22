function [L, Gt, Qt, Qt_star] = logLikelihood(model, specification, outputs, thetaD)
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

    % Inicializar la log-verosimilitud
    ll = 0;
    
    % function ll=ll_type(model,specification, d,thetaS,outputs,t)
    
    L = ll_type(model, specification, outputs, 1);
    
    for t = 2:T+1
        % Calcular Gt
        Gt(:, :, t) = calcGt(model,specification, outputs, thetaD, Gt(:, :, t-1));
        
        % Agregar un pequeño término regularizador a Gt para evitar singularidades
        reg_term = 1e-6 * eye(d);
        Gt_reg = Gt(:, :, t) + reg_term;

        % Calcular log-verosimilitud
        ll = ll_type(model, specification, outputs, t);
        L = L + ll;
        
    end
          
end
