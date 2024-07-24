LL=LogLikelihood_fgroup(model, specification, outputs, thetaD)
 % Inicialización de variables
    [i, j] = models_index(model, specification);
    rt=  outputs.returns;
    et = outputs.rotated_returns;
    T = outputs.T;
    d = outputs.d;
    thetaS=outputs.H_bar;
    P = outputs.P;
    Dt = outputs.Dt;
    Ct = outputs.Ct;
    Lambda = outputs.Lambda;
    Dt = outputs.Dt;
    Id=eye(d);

    %%%  initial inputs  %%%%
    initial_Gt=outputs.initial_Gt;
    initial_Qt_star=Id;
    delta=1;

        %fprintf('et size: %s\n', mat2str(size(et)));


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calc_all_Gts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
     Gt(:,:,1) = calcGt(model,specification,outputs, thetaD, initial_Gt)

     for t= 2:T % Adjusting loop to account for t starting from 1
        
        Gt(:, :, t) = calcGt(model,specification,outputs, thetaD, Gt(:, :, t-1)); 
            
        % Adding a small regularization term to Gt to avoid singularity
        reg_term = 1e-6 * eye(d);
        Gt_reg = Gt(:, :, t) + reg_term;
     end


        if any(isnan(Gt(:,:,t)))
            error('Gt contiene valores NaN en t=%d', t);
        end

        % fprintf('Gt size: %s\n', mat2str(size(Gt(:,:,t))));
        % disp('Gt:'); 
        % disp(Gt(:,:,t));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcQt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Qt=zeros(d,d,T)
        Qt_star=zeros(d,d,T)
        Ct=zeros(d,d,T)

     % Inicializar Qt y Qt_star si i==4
        if i == 4 % 4 is the model index for RDCC

            Qt_star(:, :, 1) = calcGt(model,specification,outputs, thetaD, initial_Qt_star); 

            for t= 2:T % Adjusting loop to account for t starting from 1
            
                Qt_star(:, :, t) = calcGt(model,specification,outputs, thetaD, Qt_star(:, :, t-1)); 
                % Adding a small regularization term to Gt to avoid singularity
                reg_term = 1e-6 * eye(d);
                Qt_star_reg = Gt(:, :, t) + reg_term;
            end   

            initial_Qt = P * sqrt(Lambda) * P' * initial_Qt_star * P * sqrt(Lambda) * P';
            initial_Ct = sqrt(initial_Qt .* Id) * initial_Qt * sqrt(inv(initial_Qt .* Id));

            % Calcular Qt y Ct para t >= 1
            for     t=1:T
                Qt(:, :, t) = P * sqrt(Lambda) * P' * Qt_star(:, :, t) * P * sqrt(Lambda) * P';
                Ct(:, :, t) = sqrt(inv(Qt(:, :, t) .* Id)) * Qt(:, :, t) * sqrt(inv(Qt(:, :, t) .* Id));
            end
        else
            Qt = [];
            Qt_star = [];
            Ct=[]
        end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ll_type_internal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fprintf('Calculando ll_type para el modelo %s en %d\n', model, t);

        function ll= ll_type_iternal(t)
          % the same as ll_type, but Ct has to be ingresed to the function as an input argument for internal use

            % Obtener la matriz U(delta)
            U =@(delta,d) delta * eye(d);
    
            switch model
                case 'RBEKK'
                    ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(t, :) * inv(Gt(:,:,t)) * et(t, :)');

                case 'OGARCH'
                    ll = -0.5 * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(t, :) * U * inv(Gt(:,:,t)) * U' * et(t, :)');

                case 'GOGARCH'
                    ll = -0.5 * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(t, :) * U * inv(Gt(:,:,t)) * U' * et(t, :)');
                case 'RDCC'
                    fprintf('Calculando RDCC, t=%d\n', t);
                    %fprintf('Dt size: %s\n', mat2str(size(Dt(:,:,t))));
                    %fprintf('Ct size: %s\n', mat2str(size(Ct(:,:,t))));

                    ll_1 = -0.5 * (d * log(2*pi) + 2*log(det(Dt(:,:,t)))) + rt(t, :) * inv(Dt(:,:,t)^2) * rt(t, :)';
                    ll_2 = -0.5 * (-et(t, :)*et(t, :)' + log(det(Ct(:,:,t))) + et(t, :)* inv(Ct(:,:,t))*et(t, :)');

                    ll = ll_1 + ll_2;
            end
        end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LogLikelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        LL = 0;

        for t= 1:T
            % Calcular log-verosimilitud
            ll = ll_type_internal(t);
            LL = LL + ll;   
        end
        end


          
     
