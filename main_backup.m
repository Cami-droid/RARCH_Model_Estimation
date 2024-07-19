% Archivo principal para ejecutar todas las funciones y generar la tabla final

% Añadir la carpeta de funciones al path
addpath('functions'); 

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');

% Extraer las columnas relevantes
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
AA = data.AA;
XOM = data.XOM;

% Calcular los log returns
log_returns_AA = diff(log(AA)) * 100;
log_returns_XOM = diff(log(XOM)) * 100;

% Combinar log returns en una matriz

log_returns = [log_returns_AA, log_returns_XOM];


% Asegurar que el vector de fechas coincida con la longitud de los log returns
dates = dates(2:end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RBEKK ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 'RBEKK'; % Cambia esto según el modelo deseado

% 1. Load and prepare data
returns = prepare_data(log_returns,model);

% 2. Rotate data
[rotated_returns, H_bar, Lambda, P] = rotate_data(returns, model);

% 3. Specificate parameters
params = specificate_params(model, 'Scalar');
params = specificate_params(model, 'Diagonal');
params = specificate_params(model, 'CP');

% 4. Conditional Covariance Estimation
cond_cov = conditional_covariance(params, model);

% 5. Generate Table
generate_table(struct('model', model, 'specification', specification, 'params', params, 'cond_cov', cond_cov));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OGARCH ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


model = 'OGARCH'; % Cambia esto según el modelo deseado

% 1. Load and prepare data
returns = prepare_data(log_returns,model);

% 2. Rotate data
[rotated_returns, H_bar, Lambda, P] = rotate_data(data, model);

% 3. optimize logLikelihood

thetaD_initial = ones(numel(Lambda_hat), 1);
 
specification='Scalar'

[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)
specification='Diagonal'
[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)
specification='CP'
[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)

% 4. Conditional Covariance Estimation
cond_cov = conditional_covariance(params, model);

% 5. Generate Table
generate_table(struct('model', model, 'specification', specification, 'params', params, 'cond_cov', cond_cov));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GOGARCH ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


model = 'GOGARCH'; % Cambia esto según el modelo deseado

% 1. Load and prepare data
returns = prepare_data(log_returns,model);

% 2. Rotate data
[rotated_returns, H_bar, Lambda, P] = rotate_data(data, model);

% 3. optimize logLikelihood

thetaD_initial = ones(numel(Lambda_hat), 1);
 
specification='Scalar'

[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)
specification='Diagonal'
[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)
specification='CP'
[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)

% 4. Conditional Covariance Estimation
cond_cov = conditional_covariance(params, model);

% 5. Generate Table
generate_table(struct('model', model, 'specification', specification, 'params', params, 'cond_cov', cond_cov));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDCC ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 'RDCC'; % Cambia esto según el modelo deseado
specification = 'Scalar'; % Cambia esto según la especificación deseada


% 1. Load and prepare data
returns = prepare_data(log_returns,model);

% 2. Rotate data
[rotated_returns, H_bar, Lambda, P] = rotate_data(data, model);

% 3. optimize logLikelihood

thetaD_initial = ones(numel(Lambda_hat), 1);
 
specification='Scalar'

[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)
specification='Diagonal'
[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)
specification='CP'
[thetaD_opt, fval, exitflag, output]=(thetaS,thetaD, thetaD_initial, rotated_returns,model,specification)

% 4. Conditional Covariance Estimation
cond_cov = conditional_covariance(params, model);

% 5. Generate Table
generate_table(struct('model', model, 'specification', specification, 'params', params, 'cond_cov', cond_cov));