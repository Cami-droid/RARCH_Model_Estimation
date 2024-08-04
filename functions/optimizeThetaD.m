function [thetaD_opt, fval, exitflag, output, L, L_marginal] = optimizeThetaD(model, specification, outputs, thetaD_initial)
%%% fval: final log-likelihood
    d=outputs.d;

    [i,j] = models_index(model, specification);
    alpha_lb=0.001;
    alpha_ub=0.999;
    beta_lb=0.85;
    beta_ub=0.99;
    lambda_cp_lb=0.001;
    lambda_cp_ub=0.99;

    %  optimization options
    options = optimset('fmincon');
    options.Display = 'off'; % 'iter'
    options.Diagnostics = 'off'; % 'on'
    options.Algorithm = 'interior-point';

    % restriccions configuration
    switch specification
        % Restricciones para el caso Scalar
        case 'Scalar'
            A = [1, 1]; % alpha + beta <= 1
            b = 1; % <=1
              lb = [0.0001, 0.9]; % alpha>=0, beta>=0
              ub = [0.99,0.99]; 
            nonlcon = [];

        % Restriccions for Diagonal specification
        case 'Diagonal'
            A = [eye(d),eye(d)];      %\alpha_{11} + \beta_{11} <= 1
                                   % \alpha_{22} + \beta_{22} <= 1
            
            b = ones(d,1);
            
            
            lb=[alpha_lb*ones(1,d),beta_lb*ones(1,d)]; % \alpha_{11},..., \alpha_{dd}, \beta_{11},..., \beta_{dd} >= 0
            ub = [alpha_ub*ones(1,d),beta_ub*ones(1,d)]; % \alpha_{11},..., \alpha_{dd}, \beta_{11},..., \beta_{dd} < 1
            nonlcon = [];

        % Restriccions for Common Persistence (CP) specification
        case 'CP'
            A = eye(d+1); % Assure that \alpha11 < 1, \alpha22 < 1, \lambda < 1
            b = ones(d+1,1);
             lb = [alpha_lb*ones(1,d), lambda_cp_lb]; % Assure that 0 < \lambda
             ub = [alpha_ub*ones(1,d), lambda_cp_ub]; % 
            nonlcon = @nonlcon;
   
    end

    if i == 3  % when the model is 'GOGARCH' add one restriction to delta
        n_rowsA = size(A, 1);
        zeros_column = zeros(n_rowsA, 1);
        A = [A, zeros_column];
       lb = [lb, 0.01];
       ub = [ub, 1];
    end

    % Definition of the negative log-likelihood function
    logLikelihoodFunc = @(thetaD) ll_engine_wrapper(model, specification, outputs, thetaD);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [thetaD_opt, fval, exitflag, output] = fmincon(logLikelihoodFunc, thetaD_initial, A, b, [], [], lb, ub, nonlcon, options);
    %fprintf('Finalizing optimizeThetaD for model %s and specification %s\n', model, specification);
    % disp(thetaD_opt);

     % Get vector L with optimized thetaD
    L = ll_engine(model, specification, outputs, thetaD_opt);

    % Getting the vector L_marginal (Txd)

    L_marginal = ll_marginal_engine(model, outputs, thetaD_opt);

    %%
    function LL = ll_engine_wrapper(model, specification, outputs, thetaD)
        L = ll_engine(model, specification, outputs, thetaD);
        LL = sum(L);
    end

    function LL_marginal=ll_marginal_engine_wrapper(model,outputs, thetaD)
        L= ll_marginal_engine(model,outputs,thetaD);
        LL= sum(LL_marginal);
    end

end
