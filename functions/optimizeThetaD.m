function [thetaD_opt, fval, exitflag, output] = optimizeThetaD(H_bar, thetaD_initial, rotated_returns, model, specification)
    % Opciones de optimización
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
    
    % Función objetivo para la log-verosimilitud
    logLikelihoodFunc = @(thetaD) logLikelihood(H_bar, thetaD, rotated_returns, model, specification);
    
    % Ejecutar la optimización
    [thetaD_opt, fval, exitflag, output] = fminunc(logLikelihoodFunc, thetaD_initial, options);
end
