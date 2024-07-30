function [thetaD_opt, fval, exitflag, output, L] = optimizeThetaD(model, specification, outputs, thetaD_initial)
%%% fval: log-verosimilitud final

    [i,j] = models_index(model, specification);
    lb=[];
    ub=[];

    % Opciones de optimización
    options = optimset('fmincon');
    options.Display = 'iter'; % 'iter'
    options.Diagnostics = 'off'; % 'on'
    options.Algorithm = 'interior-point';

    % Configuración de restricciones
    switch specification
        % Restricciones para el caso Scalar
        case 'Scalar'
            A = [1, 1]; % alpha + beta <= 1
            b = 1; % <=1
              lb = [0.0001, 0.9]; % alpha>=0, beta>=0
              ub = [0.99,0.99]; % No hay límite superior
            nonlcon = [];

        % Restricciones para el caso Diagonal
        case 'Diagonal'
            A = [
                1, 0, 1, 0;  % \alpha_{11} + \beta_{11} <= 1
                0, 1, 0, 1   % \alpha_{22} + \beta_{22} <= 1
            ];
            b = [1; 1];
         	   lb = [0.0001, 0.0001, 0.9, 0.9]; % \alpha_{11}, \alpha_{22}, \beta_{11}, \beta_{22} >= 0
               ub = [0.99,0.99,0.99,0.99];
            nonlcon = [];

        % Restricciones para el caso Common Persistence (CP)
        case 'CP'
            A = [1, 0, 0;
                 0, 1, 0;
                 0, 0, 1]; % Asegurar \alpha < 1, \beta < 1, \lambda < 1
            b = [1; 1; 1];
             lb = [0.0001, 0.0001, 0.01]; % Asegurar 0 < \lambda
             ub = [0.99,0.99,inf]; % No se necesita limite superior specific ya que \lambda < 1 está en A y b
            nonlcon = @nonlcon;
   
    end

    if i == 3  % when the model is 'GOGARCH' add one restriction to delta
        n_rowsA = size(A, 1);
        zeros_column = zeros(n_rowsA, 1);
        A = [A, zeros_column];
       lb = [lb, 0];
       ub = [ub, +1];
    end

    % Definicion de la funcion de log-verosimilitud negativa
    logLikelihoodFunc = @(thetaD) ll_engine_wrapper(model, specification, outputs, thetaD);

    % Ejecutar la optimización
    fprintf('This is thetaD_initial\n');
    disp(thetaD_initial);
    disp('This is ll_engine value at thetaD_initial\n');
    ll_engine(model, specification, outputs, thetaD_initial);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [thetaD_opt, fval, exitflag, output] = fmincon(logLikelihoodFunc, thetaD_initial, A, b, [], [], lb, ub, nonlcon, options);
    fprintf('Finalizing optimizeThetaD for model %s and specification %s\n', model, specification);
    disp(thetaD_opt);

     % Obtén el vector L con el thetaD optimizado
    L = ll_engine(model, specification, outputs, thetaD_opt);

    %%
    function LL = ll_engine_wrapper(model, specification, outputs, thetaD)
        L = ll_engine(model, specification, outputs, thetaD);
        LL = sum(L);
    end

end
