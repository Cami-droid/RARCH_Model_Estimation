% Main file to run all the functions and generate the final table
clear;clc;

% Add 'functions' folder to the path
addpath('functions'); 

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');

% Extract data from relevant columns
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
AA = data.AA; 
XOM = data.XOM;

% Calculate log returns
log_returns_AA = diff(log(AA)) * 100;
log_returns_XOM = diff(log(XOM)) * 100;

% Combine log returns in a matrix
log_returns = [log_returns_AA, log_returns_XOM];

d = size(log_returns, 2);
T = size(log_returns, 1);

% Assure that the date vector equalizes log returns' length
dates = dates(2:end);

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};
             
I = length(models);
J = length(specifications);

results(I, J) = struct('model', [], 'specification', [], 'thetaD_opt', [], 'fval', [], 'Qt', [], 'Qt_star', [],'L',[]);
outputs(I, J) = struct('model', [], 'specification', [], 'P', [], 'Lambda', [], 'H_bar', [], 'Gt', [],'returns', [], 'initials_thetaD',[],'rotated_returns', [], 'Dt', [], 'Ct', [], 'I', [], 'J', [], 'd', [], 'T', [],'L',[]);
    
initials_thetaD= { [0.02 0.91]   ,  [0.02 0.02 0.91 0.91]    ,[0.02 0.02 0.04]};
initial_delta=0.1; % only for 'GOGARCH' models
for i = 1:I
    
    for j = 1:J
        outputs(i,j).I = I;
        outputs(i,j).J = J;
        outputs(i,j).d = d;
        outputs(i,j).T = T;
        outputs(i,j).Gt = zeros(d, d, T + 1); % T+1 because the first matrix is index 0 in theory
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
        [outputs(i,j).returns, outputs(i,j).Dt] = prepare_data(model, outputs(i,j), log_returns);
        
        % Rotate data
        [outputs(i,j).rotated_returns, outputs(i,j).H_bar, outputs(i,j).Lambda, outputs(i,j).P] = rotate_data(outputs(i,j), model);

         %   disp('Unconditional Covariance Matrix H_bar valid for RBEKK, OGARCH and GOGARCH');
         %   disp(outputs(i,j).H_bar);
         %   % Matriz de eigenvectores P y matriz diagonal de eigenvalores Lambda
         %   disp('Eigenvectors matrix P:');
         %   disp(outputs(i,j).P);
         %   disp('Eigenvalues Diagnonal matrix Lambda:');
         %   disp(outputs(i,j).Lambda);

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
        
        %%%%%%%%%%%%%%%  OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [results(i,j).thetaD_opt, results(i,j).fval, exitflag, output, outputs(i,j).L] = optimizeThetaD(model, specification, outputs(i,j), initial_thetaD);
        
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

% Generate Table
generate_matlabTable;


