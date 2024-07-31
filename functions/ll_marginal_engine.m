function L_marginal = ll_marginal_engine(model, outputs, thetaD)
    % ll_marginal_engine: Calcula las log-verosimilitudes marginales para cada activo.
    % Entrada:
    %   model: String indicando el modelo (e.g., 'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC')
    %   outputs: Estructura con datos de entrada (r_t, etc.)
    %   thetaD: Vector de parámetros dinámicos
    % Salida:
    %   L_marginal: Matriz de tamaño Txd con las log-verosimilitudes marginales

    % Variables iniciales
    rt = outputs.returns; % Retornos originales, matriz Txd
    et = outputs.rotated_returns; % Retornos rotados, matriz Txd
    T = outputs.T; % Número de observaciones
    d = outputs.d; % Número de activos
    Dt = outputs.Dt; % Matrices diagonales con desviaciones estándar, matriz dxdxT
    Lambda = outputs.Lambda; % Valores propios para OGARCH y GOGARCH
    P = outputs.P; % Matriz de autovectores para OGARCH y GOGARCH
    delta = thetaD(end); % Parámetro de delta para OGARCH y GOGARCH

    % Inicializar matrices de salida
    L_marginal = zeros(T, d);

    % Calcular las log-verosimilitudes marginales
    switch model
        case 'RBEKK'
            for i = 1:d
                for t = 1:T
                    sigma2_t = outputs.Gt(i, i, t); % Varianza condicional
                    L_marginal(t, i) = -0.5 * (log(2 * pi) + log(sigma2_t) + (et(t, i)^2) / sigma2_t);
                end
            end

        case {'OGARCH', 'GOGARCH'}
            U = @(delta, d) delta * eye(d);
            for i = 1:d
                for t = 1:T
                    sigma2_t = U(delta, d) * outputs.Gt(:,:,t) * U(delta, d)';
                    L_marginal(t, i) = -0.5 * (log(2 * pi) + log(sigma2_t(i, i)) + (et(t, i)^2) / sigma2_t(i, i));
                end
            end

        case 'RDCC'
            for i = 1:d
                for t = 1:T
                    sigma2_t = Dt(i, i, t)^2; % Varianza condicional diagonal
                    L_marginal(t, i) = -0.5 * (log(2 * pi) + log(sigma2_t) + (rt(t, i)^2) / sigma2_t);
                end
            end
    end
end 