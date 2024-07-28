% Main file to run all the functions and generate the final table
clear;clc;

% Add 'functions' folder to the path
addpath('functions'); 

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');
height(data)
% Extract data from relevant columns, Calculate log returns and Combine them in a matrix
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
T=height(data)-1
log_retrurns=zeros(1,28)
log_ret_AA=diff(log(data.AA))*100;
log_ret_AXP=diff(log(data.AXP))*100;
log_ret_BA=diff(log(data.BA))*100;
log_ret_BAC=diff(log(data.BAC))*100;
log_ret_CAT=diff(log(data.CAT))*100;
log_ret_CSCO=diff(log(data.CSCO))*100;
log_ret_CVX=diff(log(data.CVX))*100;
log_ret_DD=diff(log(data.DD))*100;
log_ret_DIS=diff(log(data.DIS))*100;
log_ret_GE=diff(log(data.GE))*100;
log_ret_HD=diff(log(data.HD))*100;
log_ret_HPQ=diff(log(data.HPQ))*100;
log_ret_IBM=diff(log(data.IBM))*100;
log_ret_INTC=diff(log(data.INTC))*100;
log_ret_JNJ=diff(log(data.JNJ))*100;
log_ret_JPM=diff(log(data.JPM))*100;
log_ret_KO=diff(log(data.KO))*100;
log_ret_MCD=diff(log(data.MCD))*100;
log_ret_MMM=diff(log(data.MMM))*100;
log_ret_MRK=diff(log(data.MRK))*100;
log_ret_MSFT=diff(log(data.MSFT))*100;
log_ret_PFE=diff(log(data.PFE))*100;
log_ret_PG=diff(log(data.PG))*100;
log_ret_T=diff(log(data.T))*100;
log_ret_TRV=diff(log(data.TRV))*100;
log_ret_VZ=diff(log(data.VZ))*100;
log_ret_WMT=diff(log(data.WMT))*100;
log_ret_XOM=diff(log(data.XOM))*100;

% Combine log returns in a matrix
log_returns = [log_ret_AA, log_ret_XOM, log_ret_BA, log_ret_BAC,log_ret_CAT log_ret_CSCO log_ret_CVX log_ret_DD log_ret_DIS log_ret_GE log_ret_HD log_ret_HPQ log_ret_IBM log_ret_INTC log_ret_JNJ log_ret_JPM log_ret_KO log_ret_MCD log_ret_MMM log_ret_MRK log_ret_MSFT log_ret_PFE log_ret_PG log_ret_T log_ret_TRV log_ret_VZ log_ret_WMT log_ret_XOM ] ;
d = size(log_returns, 2);
T = size(log_returns, 1);

% Assure that the date vector equalizes log returns' length
dates = dates(2:end);

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};
             
I = length(models);
J = length(specifications);

results(I, J) = struct('model', [], 'specification', [], 'thetaD_opt', [], 'fval', [], 'Qt', [], 'Qt_star', [],'L',[]);
outputs(I, J) = struct('model', [], 'specification', [], 'P', [], 'Lambda', [], 'H_bar', [], 'Gt', [],'returns', [], 'initials_thetaD',[],'rotated_returns', [], 'Dt', [], 'Ct', [], 'I', [], 'J', [], 'd', [], 'T', []);
 
init_alpha_scalar=0.05
init_beta_scalar=0.9
init_alpha_diag=ones(1,d)*init_alpha_scalar
init_beta_diag=ones(1,d)*init_beta_scalar


for i = 1:I
    for j = 1:J
        outputs(i,j).I = I;
        outputs(i,j).J = J;
        outputs(i,j).d = d;
        outputs(i,j).T = T;
        outputs(i,j).Gt = zeros(d, d, T + 1); % T+1 because the first matrix is index 0 in theory
        outputs(i,j).initial_Gt = eye(d);


        if i==3
        % Only for GOGARCH
        outputs(3,j).initial_delta=1;
        
        outputs(i,j).initials_thetaD = { [init_alpha_scalar init_beta_scalar outputs(3,j).initial_delta]      , [init_alpha_diag init_beta_diag outputs(3,j).initial_delta]       , [init_alpha_diag 0.993 outputs(3,j).initial_delta]};
        
        else
        outputs(i,j).initials_thetaD = { [init_alpha_scalar init_beta_scalar]      , [init_alpha_diag init_beta_diag]       , [init_alpha_diag 0.3]};
        end
        
    end
end
%% 

for i = 1:I
    model = models{i};
    fprintf('Estimating model: %s\n', model);
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
        initial_thetaD = outputs(i,j).initials_thetaD{j};

        fprintf('The initials thetaDs are: %s\n', mat2str(outputs(i,j).initials_thetaD{j}));
        
        % Estimate thetaD parameters
        outputs(i,j).Gt= calc_all_Gts(model, specification, outputs(i,j), outputs(i,j).initials_thetaD{j}, outputs(i,j).initial_Gt);
        
        [results(i,j).thetaD_opt, results(i,j).fval, exitflag, output] = optimizeThetaD(model, specification, outputs(i,j), initial_thetaD);
        
        fprintf('The optimal thetaDs found are: %s\n', mat2str(results(i,j).thetaD_opt));
        disp(results(i,j).thetaD_opt)
        fprintf('LogLikelihood value: %s\n', mat2str(results(i,j).fval));

        % calculate Gt at the optimum thetaD_opt

        Id=eye(d);
        output(i,j).Gt=calcGt(model, specification, outputs(i,j),results(i,j).thetaD_opt,Id);
        
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

% Generar la tabla y guardarla en un archivo .txt
thetaD_table_matlab = generate_pdfTable(results);

% Crear y compilar el archivo LaTeX

% Define the paths
results_dir = 'D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\results';
report_file = fullfile(results_dir, 'report.tex');

% Compile the LaTeX file to PDF, specifying the output directory
command = sprintf('pdflatex -output-directory="%s" "%s"', results_dir, report_file);
status = system(command);

% Check if the command was successful
if status == 0
    disp('PDF generated successfully.');
else
    disp('Error generating PDF.');
end 