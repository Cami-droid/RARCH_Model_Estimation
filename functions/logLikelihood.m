function [LL, Gt, Qt, Qt_star L_matrix] = logLikelihood(model, specification, outputs, thetaD)
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
    
    
    LL = 0;

    for t= 1:T
        % Calcular log-verosimilitud
        ll = ll_type(model, specification, outputs, t);
        LL = LL + ll;        
    end
          
end
