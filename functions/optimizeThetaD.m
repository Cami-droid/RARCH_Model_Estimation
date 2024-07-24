function [thetaD_opt, fval, exitflag, output, L] = optimizeThetaD(model, specification, outputs, thetaD_initial)
    models_index(model, specification);

    % Opciones de optimización
    options = optimset('fmincon');
    options.Display = 'off'; % 'iter'
    options.Diagnostics = 'off'; % 'on'
    options.Algorithm = 'interior-point';

    % Configuración de restricciones
    switch specification
        % Restricciones para el caso Scalar
        case 'Scalar'
            A = [1, 1]; % alpha + beta <= 1
            b = 1; % <=1
            lb = [0, 0]; % alpha>=0, beta>=0
            ub = []; % No hay límite superior
            nonlcon = [];

        % Restricciones para el caso Diagonal
        case 'Diagonal'
            A = [
                1, 0, 1, 0;  % \alpha_{11} + \beta_{11} <= 1
                0, 1, 0, 1   % \alpha_{22} + \beta_{22} <= 1
            ];
            b = [1; 1];
            lb = [0, 0, 0, 0]; % \alpha_{11}, \alpha_{22}, \beta_{11}, \beta_{22} >= 0
            ub = [];
            nonlcon = [];

        % Restricciones para el caso Common Persistence (CP)
        case 'CP'
            A = [1, 0, 0;
                 0, 1, 0;
                 0, 0, 1]; % Asegurar \alpha < 1, \beta < 1, \lambda < 1
            b = [1; 1; 1];
            lb = [0, 0, 0]; % Asegurar 0 < \lambda
            ub = []; % No se necesita límite superior específico ya que \lambda < 1 está en A y b
            nonlcon = @nonlcon;

    end

    % Definición de la función de log-verosimilitud negativa
    logLikelihoodFunc = @(thetaD) -ll_engine(model, specification, outputs, thetaD);

    % Ejecutar la optimización
    [thetaD_opt, fval, exitflag, output] = fmincon(logLikelihoodFunc, thetaD_initial, A, b, [], [], lb, ub, nonlcon, options);

    disp(thetaD_opt);

    % Log-Verosimilitud final
    L = -fval; % Recordar que fval es la log-verosimilitud negativa, por lo que invertimos el signo
end
