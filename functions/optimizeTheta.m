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
                [theta_vec, fval, Gt, VCV, Scores] = rarch(et, 1, 1, specification, '2-stage');

            case {'OGARCH', 'GOGARCH'}
                [theta_vec, fval, Gt, VCV, Scores] = gogarch(et, 1, 1, [], model);

            case 'RDCC'
                [theta_vec, fval, Gt, VCV, Scores] = rdcc(rt, [], 1, [], 1, 1, 0, 1, 2, '2-stage', 'Scalar', [], []);

            otherwise
                error('Unknown model: %s', model);
        end
    catch ME
        fprintf('Error optimizing model %s with specification %s: %s\n', model, specification, ME.message);
    end

    % Additional checks or calculations can be added here if needed

end

