% Main file to run all the functions and generate the final table
clear;clc;
global outputs globals_var results


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

% Combine log returns in a matriz
log_returns = [log_returns_AA, log_returns_XOM];


d=size(log_returns,2);
T=size(log_returns,1);


% Assure  that the date vector equalize log returns's length

dates = dates(2:end);

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};
initials_thetaD = {[0.2 0.05], [0.05 0.05 0.2 0.2],[0.05 0.2 0.25]};

I=length(models);
J=length(specifications);

globals_var=struct('I',I,'J',J,'d',[],'T',[])
results(I, J) = struct('model', [], 'specification', [], 'params', []);%, 'M_params', [], 'std_errors', [], 'totalLL', [], 'Copula', [], 'Margin', []);
outputs(I,J)=struct('model',[],'specification',[],'P',[],'Lambda',[],'H_bar',[],'Gt',[],'returns',[],'rotated_returns',[],'Dt',[],'Ct',[])

for i = 1:length(models)
    model = models{i};
    fprintf('Estimating model: %s\n', model);
    % Load and prepare data
    [returns,Dt] = prepare_data(log_returns, model);
        
    % Rotate data
    [rotated_returns, H_bar, Lambda, P] = rotate_data(returns, model);
    
    for j = 1:length(specifications)
        specification = specifications{j};
        fprintf('Estimating model: %s\n', model);
        fprintf('Parameter specification: %s\n', specification);
        
        % Optimize logLikelihood
        thetaD_initial= initials_thetaD{j};
        fprintf('The initials thetaD are: %s\n', mat2str(initials_thetaD{j}));
        
        [thetaD_opt, fval, exitflag, output] = optimizeThetaD(H_bar, thetaD_initial, rotated_returns, model, specification);
                      
                   
        % Generate Table
        
        % Create the global structures
        results(i,j) = struct('model', model, 'specification', specification, 'params', thetaD_opt,'totalLL',fval);
        outputs(i,j)=struct('model',model,'specification',specification,'P',P,'Lambda',Lambda,'H_bar',H_bar,'Gt',Gt,'returns',returns,'rotated_returns',rotated_returns)

    % % Call the function
    %  generate_table(results(i,j));
         
    end
end
