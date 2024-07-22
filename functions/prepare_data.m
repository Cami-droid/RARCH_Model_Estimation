function  [returns, Dt]= prepare_data(log_returns,model)
 % compute standarized returns series and diagonal variance matrices at every t

    global d;
    global T;
    
    % Calcular la media de los log returns
    mean_log_returns = mean(log_returns);

    % Calcular los residuales (media cero)
    demeaned_returns = log_returns - mean_log_returns;
    
    % Preparar los datos segun el modelo
    switch model
        case 'RBEKK'
            % Specific preparation for RBEKK.
            prepared_data=demeaned_returns;
            
                        
        case 'OGARCH'
            % Specific preparation for OGARCH

            prepared_data=demeaned_returns;
                        
        case 'GOGARCH'
            % Specific preparation for GOGARCH

            prepared_data=demeaned_returns;
                        
        case 'RDCC'
            % Specific preparation for RDCC

            std_returns = zeros(size(demeaned_returns));
            
            global Dt
            Dt = zeros(d,d,T);
            
            cond_var = zeros(size(demeaned_returns));
                
            for i = 1:d
                garch_model = garch('GARCHLags', 1, 'ARCHLags', 1);
                garch_fit = estimate(garch_model, demeaned_returns(:, i), 'display', 'off');
                cond_var(:, i) = infer(garch_fit, demeaned_returns(:, i));
                std_returns(:, i) = demeaned_returns(:, i) ./ sqrt(cond_var(:, i));
                summarize(garch_fit);

            end
                
            for t=1:T
                    Dt(:,:,t)=diag(cond_var(t,:));
            end
            
            prepared_data=std_returns;
                         
                        
    end
    
    % Devolver los datos preparados
    returns = prepared_data;
    
end 