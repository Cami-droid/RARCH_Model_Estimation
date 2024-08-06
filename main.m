% Iniciar registro de errores y salida en un archivo de texto
diary('logfile.txt');
diary on;

try
    % Main file to run all the functions and generate the final table

    % Add 'functions' folder to the path
    addpath('functions'); 

    models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
    specifications = {'Scalar', 'Diagonal', 'CP'};
             
    I = length(models);
    J = length(specifications);

    results(I, J) = struct('model', [], 'specification', [], 'thetaM',[],'thetaD_opt', [], 'fval', [], 'Qt', [], 'Qt_star', [],'L',[],'LL_marginal',[],'LL_copula',[]);
    outputs(I, J) = struct('model', [], 'specification', [], 'P', [], 'Lambda', [], 'H_bar', [], 'Gt', [],'returns', [],'std_returns',[], 'initials_thetaD',[],'rotated_returns', [], 'Dt', [], 'Ct', [], 'I', [], 'J', [], 'd', [], 'T', [],'L',[],'LL_marginal',[],'LL_copula',[]);
    alpha_init=0.01;
    beta_init=0.98;


    initials_thetaD= { [alpha_init, beta_init], [alpha_init*ones(1,d),beta_init*ones(1,d)],[alpha_init*ones(1,d), alpha_init+beta_init]};

    initial_delta=0.1; % only for 'GOGARCH' models

    
    for i = 1:I

        for j = 1:J
            outputs(i,j).I = I;
            outputs(i,j).J = J;
            outputs(i,j).d = d;
            outputs(i,j).T = T;
            outputs(i,j).Gt = zeros(d, d, T); % T+1 because the first matrix is index 0 in theory
            outputs(i,j).initial_Gt = eye(d);
            outputs(i,j).initials_thetaD=initials_thetaD{j};

            if i==3
                outputs(i,j).initials_thetaD=[initials_thetaD{j} initial_delta];
            end

        end
    end
    clear delta

    %% 

    for i = 1:I
        model = models{i};
        fprintf('*********************************** Estimating model: %s ****************************************\n', model);
        for j = 1:J
            % Load and prepare data
            % Calculate the mean of  log returns
            mean_log_returns = mean(log_returns);

            % Calculate residuals rt (zero mean)
            outputs(i,j).returns= log_returns - mean_log_returns;

            [outputs(i,j).std_returns, outputs(i,j).Dt, results(i,j).thetaM] = prepare_data(model, outputs(i,j), log_returns);

            % Rotate data
            [outputs(i,j).rotated_returns, outputs(i,j).H_bar, outputs(i,j).Lambda, outputs(i,j).P,  outputs(i,j).PI_bar, outputs(i,j).Lambda_C, outputs(i,j).P_C] = rotate_data(outputs(i,j), model);

        end

        for j = 1:J
            specification = specifications{j};
            outputs(i,j).model = model;
            outputs(i,j).specification = specification;

            fprintf('Estimating model: %s\n', model);
            fprintf('Parameter specification: %s\n', specification);

            % Optimize logLikelihood
            initial_thetaD = outputs(i,j).initials_thetaD;

            fprintf('The initials thetaDs are: %s\n', mat2str(outputs(i,j).initials_thetaD));

            % Estimate thetaD parameters
            outputs(i,j).Gt= calc_all_Gts(model, specification, outputs(i,j), outputs(i,j).initials_thetaD, outputs(i,j).initial_Gt);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [results(i,j).thetaD_opt, results(i,j).fval, exitflag, output, outputs(i,j).L,outputs(i,j).L_marginal] = optimizeTheta(model, specification, outputs(i,j), initial_thetaD);

            % Storing LL_marginals in results

            %results(i,j).LL_marginal=sum(outputs(i,j).L_marginal)';
            %results(i,j).LL_copula=results(i,j).fval-sum(results(i,j).LL_marginal);


            fprintf('The optimal thetaDs found are: %s\n', mat2str(results(i,j).thetaD_opt));
            fprintf('LogLikelihood value: %s\n', mat2str(results(i,j).fval));

            % calculate Gt at the optimum thetaD_opt

            Id=eye(d);
            output(i,j).Gt=calc_all_Gts(model, specification, outputs(i,j),results(i,j).thetaD_opt,Id);

            % Calculate Qt, Qt_star and Ct

            [outputs(i,j).Qt, outputs(i,j).Qt_star outputs(i,j).Ct] = calcQt(model, specification, outputs(i,j), results(i,j).thetaD_opt);

            % Store the results
            results(i,j).model = model;
            results(i,j).specification = specification;
        end
    end
    %% Table Generation
    generate_matlabTable;

catch ME
    % Guardar información del error
    fprintf('An error occurred: %s\n', ME.message);
    fprintf('Error occurred in:\n');
    for k = 1:length(ME.stack)
        fprintf('File: %s, Line: %d\n', ME.stack(k).file, ME.stack(k).line);
    end
end

% Terminar registro de errores y salida
diary off;
