function [theta_vec, fval, Gt, VCV, Scores] = optimizeTheta(model, specification, outputs, thetaD_initial)
    % Function to optimize theta parameters for different models
    % Inputs:
    %   model          - The model name ('RBEKK', 'OGARCH', 'GOGARCH', 'RDCC')
    %   specification  - The specification type
    %   outputs        - The structure containing necessary data and parameters
    %   thetaD_initial - Initial thetaD parameters
    % Outputs:
    %   theta_vec      - Optimized theta parameters
    %   fval           - Final log-likelihood value
    %   Gt             - Estimated Gt matrices
    %   VCV            - Variance-covariance matrix
    %   Scores         - Scores

    % Extract necessary data
    d = outputs.d;
    et = outputs.rotated_returns;
    rt = outputs.returns;

    %  optimization options
    options = optimset('fmincon');
    options.Display = 'off'; % 'iter'
    options.Diagnostics = 'off'; % 'on'
    options.Algorithm ='interior-point';%'sqp'

    % Initialize outputs
    
    theta_vec = [];
    fval = [];
    Gt = [];
    VCV = [];
    Scores = [];

    % Switch case to call the appropriate optimization function based on the model
    try
        switch model
            case 'RBEKK'
                [theta_vec, fval, Gt, VCV, Scores] = rarch(et, 1, 1, specification, '2-stage',[],options);
                %rarch(data,p,q,type,method,startingVals,options)

            case {'OGARCH', 'GOGARCH'}
                [theta_vec, fval, Gt, VCV, Scores] = gogarch(et, 1, 1, [], model,[],options);
                %gogarch(data,p,q,gjrType,type,startingVals,options)

            case 'RDCC'
                [theta_vec, fval, Gt, VCV, Scores]             = dcc(rt, []     , 1,[],1,1,0,1,2,'2-stage', 'None', [], options);
                %[parameters, ll ,Ht, VCV, scores, diagnostics]=dcc(data,dataAsym,m,l,n,p,o,q,gjrType,method,composite,startingVals,options)

            otherwise
                error('Unknown model: %s', model);
        end
    catch ME
        fprintf('Error optimizing model %s with specification %s: %s\n', model, specification, ME.message);
    end

    % Additional checks or calculations can be added here if needed

end 