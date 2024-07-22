function [L, Gt, Qt, Qt_star] = logLikelihood(model, specification, outputs, thetaD)
    % Inicialización de variables
    [i, j] = models_index(model, specification);
    rotated_returns = outputs.rotated_returns;
    T = outputs.T;
    d = outputs.d;
    thetaS=outputs.H_bar;
    P = outputs.P;
    Gt =outputs.Passenger_Gt;
    Lambda = outputs.Lambda;
    Dt = outputs.Dt;

    % Inicializar la log-verosimilitud
    ll = 0;
    
    % function ll=ll_type(model,specification, d,thetaS,outputs,t)
    
    L = ll_type(model, specification, outputs, 1);
    
    for t_count = 2:T+1
        % Calcular log-verosimilitud
        ll = ll_type(model, specification, outputs, t_count-1);
        L = L + ll;
        
    end
          
end
