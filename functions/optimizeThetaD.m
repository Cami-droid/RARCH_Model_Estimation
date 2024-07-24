function [thetaD_opt, fval, exitflag, output, L] = optimizeThetaD(model,specification, outputs, thetaD_initial)

    models_index(model, specification)
    %  fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    
    % input arguments:
     
    %       thetaS is the H_bar estimator, i.e, the Static estimator
    %       thetaD_initial is the arbitrary initial vector for the initial dynamic parameter estimator vector
    %       rotated_returns 
    %       model could be RBEKK, OGARCH, GOGARCH or RDCC
    %       specification can be 'Scalar','Diagonal','CP'
        
    % Optimization settings
    
    options = optimset('fmincon');
    options.Display = 'off'; % 'iter'
    options.Diagnostics = 'off'; % 'on'
    options.Algorithm = 'interior-point';
    %%   
    %     Restrictions
    % 
    %     (1-lambda_cp)Id is positive definite.
    %     max lambda_cp_ij < 1, for all i, j = 1, . . . , d.
    %     lambda_cp_ii = alpha_ii + beta_ii < 1
    % 
   
    switch specification
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scalar restrictions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'Scalar'
            A = [1, 1]; % alpha + beta <= 1
            b = 1; % <=1
            lb = [0, 0]; % alpha>=0, beta>=0
            ub = []; % There is no upper bound
            nonlcon=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Diagonal restrictions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 'Diagonal'
            A = [
            1, 0, 1, 0;  % \alpha_{11} + \beta_{11} <= 1
            0, 1, 0, 1   % \alpha_{22} + \beta_{22} <= 1
            ];
            b = [1; 1];
            lb = [0, 0, 0, 0]; % \alpha_{11}, \alpha_{22}, \beta_{11}, \beta_{22} >= 0
            ub = [];
            nonlcon=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Common Persistence restrictions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'CP'
              A = [1, 0, 0;
                  0, 1, 0;
                  0, 0, 1]; % To assure \alpha < 1, \beta < 1, \lambda < 1

              b = [1; 1; 1];

              lb = [0, 0, 0]; % To assure 0 < \lambda

              ub = []; % We don't need an specific upper bound as \lambda < 1 is driven in A y b
              nonlcon=@nonlcon;

    end
    
    % Objective function for optimization

    logLikelihoodFunc = @(thetaD) LogLikelihood_group(model, specification, outputs, thetaD);
    
    % Run the optimization  x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
                           
    [thetaD_opt, fval, exitflag, output] = fmincon(logLikelihoodFunc, thetaD_initial, A, b, [], [], lb, ub, nonlcon, options);
    
    disp(thetaD_opt)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%results.(i,j).L=ll;
end


