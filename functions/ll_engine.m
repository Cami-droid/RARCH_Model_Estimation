function L = ll_engine(model, specification, outputs, thetaD)
    % devuelve un vector columna L con las log-verosimilitudes en cada periodo t=1...T
    % recuerda que en la entrada de outputs en realidad entra outputs(i,j)
    % Inicialización de variables
    % Imprime la frase y el vector en un solo renglón
    %sentence='trying ';
    %fprintf('%s ', sentence);
    %fprintf('%d ', thetaD);
    %fprintf('\n'); % line jump


    [i, j] = models_index(model, specification);
    rt = outputs.returns;
    et = outputs.rotated_returns;
    T = outputs.T;
    d = outputs.d;
    thetaS = outputs.H_bar;
    P = outputs.P;
    Dt = outputs.Dt;
    Ct = outputs.Ct;
    Lambda = outputs.Lambda;
    Id = eye(d);

    % Inicializar inputs iniciales

    initial_Gt = outputs.initial_Gt;
    initial_Qt_star = Id;
    delta = thetaD(end); %%%%%%%%%%%%%%%%%%%% falta desarrollar

    % Inicializar matrices de salida
    Gt = zeros(d, d, T);
    Gt(:,:,1) = calcGt(model, specification, outputs, thetaD, initial_Gt,1);
 
    % Calculate Gt for t >= 2
    for t = 2:T
        Gt(:,:,t) = calcGt(model, specification, outputs, thetaD, Gt(:,:,t-1),t);
        
        % Agregar un término de regularización pequeño a Gt para evitar singularidad
        reg_term = 1e-6 * eye(d);
        Gt(:,:,t) = Gt(:,:,t) + reg_term;
    end

    % Verificar si hay valores NaN en Gt
    if any(isnan(Gt), 'all')
        error('Gt contiene valores NaN');
    end

    % Inicializar Qt y Qt_star si el modelo es RDCC
    if i == 4 % 4 es el índice del modelo para RDCC
        Qt_star = zeros(d, d, T);
        Qt_star(:,:,1) = initial_Qt_star;

        for t = 2:T
            Qt_star(:,:,t) = calcGt(model, specification, outputs, thetaD, Qt_star(:,:,t-1),t);
            Qt_star(:,:,t) = Qt_star(:,:,t) + reg_term;
        end

        initial_Qt = P * sqrt(Lambda) * P' * initial_Qt_star * P * sqrt(Lambda) * P';
        initial_Ct = sqrt(initial_Qt .* Id) * initial_Qt * sqrt(inv(initial_Qt .* Id));

        Qt = zeros(d, d, T);
        Ct = zeros(d, d, T);
        
        for t = 1:T
            Qt(:,:,t) = P * sqrt(Lambda) * P' * Qt_star(:,:,t) * P * sqrt(Lambda) * P';
            Ct(:,:,t) = sqrt(inv(Qt(:,:,t) .* Id)) * Qt(:,:,t) * sqrt(inv(Qt(:,:,t) .* Id));
        end
    else
        Qt = [];
        Qt_star = [];
        Ct = [];
    end

        % Calcular el vector columna de log-verosimilitud en cada período t
        L=zeros(T,1);
     for t = 1:T
        ll = ll_type_internal(t);
        L(t,1)=ll;
    end


% Definir la función ll_type_internal
    function ll = ll_type_internal(t)
        U =@(delta,d) delta * eye(d);
        
        switch model
            case 'RBEKK'
                ll = -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(t, :) * inv(Gt(:,:,t)) * et(t, :)');
            case 'OGARCH'
                ll = -0.5 * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(t, :) * U(delta, d) * inv(Gt(:,:,t)) * U(delta, d)' * et(t, :)');
            case 'GOGARCH'
                ll = -0.5 * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(t, :) * U(delta, d) * inv(Gt(:,:,t)) * U(delta, d)' * et(t, :)');
            case 'RDCC'
                ll_1 = -0.5 * (d * log(2*pi) + 2*log(det(Dt(:,:,t)))) + rt(t, :) * inv(Dt(:,:,t)^2) * rt(t, :)';
                ll_2 = -0.5 * (-et(t, :)*et(t, :)' + log(det(Ct(:,:,t))) + et(t, :) * inv(Ct(:,:,t)) * et(t, :)');
                ll = ll_1 + ll_2;
        end
    end

end
