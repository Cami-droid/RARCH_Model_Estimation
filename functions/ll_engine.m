function L = ll_engine(model, specification, outputs, thetaD)
    % returns the column vector L with log-likelihoods en each period t=1...T
    % remind that in the input argument for 'outputs', actually is outputs(i,j)
    
    [i, j] = models_index(model, specification);
    rt = outputs.returns'; % transposed to fit the paper formula
    Et= outputs.std_returns'; % transposed to fit the paper formula
    et = outputs.rotated_returns';% transposed to fit the paper formula
    T = outputs.T;
    d = outputs.d;
    thetaS = outputs.H_bar;
    PI_bar=outputs.PI_bar;
    P_C=outputs.P_C;
    Lambda_C=outputs.Lambda_C;
    P = outputs.P;
    Dt = outputs.Dt;
    Ct = outputs.Ct;
    Lambda = outputs.Lambda;
    Id = eye(d);



    % Add a small regulation term to Gt to avoid singularity
    reg_term = 1e-6 * eye(d);

    % Inicialize initial inputs

    initial_Gt = outputs.initial_Gt;
    initial_Qt_star = Id;
    delta = thetaD(end); %%%%%%%%%%%%%%%%%%%% falta desarrollar

    % Inicialize output matrices
    Gt = zeros(d, d, T);
    Gt(:,:,1) = calcGt(model, specification, outputs, thetaD, initial_Gt,1);
 
    % Calculate Gt for t >= 2
    for t = 2:T
        Gt(:,:,t) = calcGt(model, specification, outputs, thetaD, Gt(:,:,t-1),t);             
    end

    % Verify NaN values in Gt
    if any(isnan(Gt), 'all')
         error('Gt contains NaN values');
    end

    % Inicialize Qt y Qt_star when the model is RDCC

    Lambda_C_inv_sqrt = diag(1 ./ sqrt(max(diag(Lambda_C))));
                            
    % specific rotation for RDCC
     PI_bar_inv_sqrt = P_C * Lambda_C_inv_sqrt * P_C';        


    if i == 4 % 4 is the model index for RDCC
        Qt_star = zeros(d, d, T);
        Qt_star(:,:,1) = initial_Qt_star;

        for t = 2:T
            Qt_star(:,:,t) = calcGt(model, specification, outputs, thetaD, Qt_star(:,:,t-1),t) + reg_term;           
        end

        initial_Qt = P_C * sqrt(Lambda_C) * P_C' * initial_Qt_star * P_C * sqrt(Lambda_C) * P_C';
        initial_Ct = sqrt(initial_Qt .* Id) * initial_Qt * sqrt(inv(initial_Qt .* Id));

        Qt = zeros(d, d, T);
        Ct = zeros(d, d, T);
        
        for t = 1:T
            Qt(:,:,t) = P_C * sqrt(Lambda_C) * P_C' * Qt_star(:,:,t) * P_C * sqrt(Lambda_C) * P_C';
            Ct(:,:,t) = sqrt(inv(Qt(:,:,t) .* Id)) * Qt(:,:,t) * sqrt(inv(Qt(:,:,t) .* Id));
        end
    else
        Qt = [];
        Qt_star = [];
        Ct = [];
    end

        % Calculate the log-likelihood column vector x in every t period 
        L=zeros(T,1);
     for t = 1:T
        ll = ll_type_internal(t);
        L(t,1)=ll;
    end

% Define the ll_type_internal function
    function ll = ll_type_internal(t)
        U =@(delta,d) delta * eye(d);
        
        switch model
            case 'RBEKK'
                ll = (-1/2) * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(:,t)'* inv(Gt(:,:,t))* et(:,t));
                
            case 'OGARCH'
                ll = (-1/2) * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(:,t)'*inv(Gt(:,:,t))*et(:,t));
                
            case 'GOGARCH'
                ll = (-1/2) * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(:,t)' * U(delta, d) * inv(Gt(:,:,t)) * U(delta, d)' * et(:,t));
               
            case 'RDCC'
                ll_M = (-1/2) * (d * log(2*pi) + 2*log(det(Dt(:,:,t))) + rt(:, t)'* (Dt(:,:,t)^(-2))*rt(:, t));
                ll_C = (-1/2) * (-Et(:,t)'*Et(:,t) + log(det(Ct(:,:,t))) + Et(:,t)'* (Ct(:,:,t)^(-1))*Et(:,t));

                ll =  ll_M+ll_C;
        end
    end

end
