% Archivo principal para ejecutar todas las funciones y generar la tabla final

% Añadir la carpeta de funciones al path
addpath('functions'); 

% Cargar los datos
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

% Modelo y especificaciones a evaluar
modelos = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
especificaciones = {'Scalar', 'Diagonal', 'CP'};

% Iterar sobre cada modelo y especificación
for m = 1:length(modelos)
    model = modelos{m};
    for s = 1:length(especificaciones)
        specification = especificaciones{s};
        
        % Paso 1: Preparar datos
        returns = prepare_data(log_returns, model);

        % Paso 2: Rotar datos
        [rotated_returns, H_bar, Lambda, P] = rotate_data(returns, model);

        % Paso 3: Optimizar logLikelihood
        thetaD_initial = ones(numel(Lambda), 1);
        [thetaD_opt, fval, exitflag, output] = optimizeThetaD(H_bar, thetaD_initial, rotated_returns, model, specification);

        % Paso 4: Estimación de la Covarianza Condicional
        cond_cov = conditional_covariance(thetaD_opt, rotated_returns, H_bar, model, specification);

        % Paso 5: Generar Tabla
        generate_table(struct('model', model, 'specification', specification, 'params', thetaD_opt, 'cond_cov', cond_cov));
    end
end
