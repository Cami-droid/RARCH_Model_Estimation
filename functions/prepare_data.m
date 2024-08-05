function  [std_returns,Dt, thetaM]= prepare_data(model,outputs,log_returns)
 % compute standarized returns series and diagonal variance matrices at every t
 d=outputs.d;
 T=outputs.T;
 rt=outputs.returns;
    

    thetaM=NaN(d*2,1);  
    % Prepare data according to the model
    switch model
        case 'RBEKK'
            % Specific preparation for RBEKK.
            prepared_data=[]; % this model doesn't work with standarized returns
            Dt=NaN;
            thetaM=[];
                                    
        case 'OGARCH'
            % Specific preparation for OGARCH

            prepared_data=[]; % this model doesn't work with standarized returns
            Dt=NaN;
            thetaM=[];
                        
        case 'GOGARCH'
            % Specific preparation for GOGARCH
            %for i=1:d
           %    garch_model = garch('GARCHLags', 1, 'ARCHLags', 1);
           %     garch_fit = estimate(garch_model, rt(:, i), 'display', 'off');
           %     cond_var(:, i) = infer(garch_fit, rt(:, i));
           %     std_returns(:, i) = rt(:, i) ./ sqrt(cond_var(:, i));
           %  end
           % prepared_data=std_returns;
           prepared_data=[]; % this model doesn't work with standarized returns
            Dt=NaN;
            thetaM=[];
                  
        case 'RDCC'
            % Specific preparation for RDCC

            std_returns = zeros(size(rt));
            Dt = zeros(d,d,T);
            
            cond_var = zeros(size(rt));
            thetaM=zeros(d*2,1); 
            for i = 1:d

                omega = 0.01; % Valor inicial para omega
                alpha = 0.1; % Valor inicial para alpha
                beta = 0.8;  % Valor inicial para beta
                garch_model = garch('GARCHLags', 1, 'ARCHLags', 1);
                garch_fit = estimate(garch_model, rt(:, i), 'display', 'off');
                cond_var(:, i) = infer(garch_fit, rt(:, i));
                std_returns(:, i) = rt(:, i) ./ sqrt(cond_var(:, i));
                %summarize(garch_fit);

                omega_est = garch_fit.Constant;
                alpha_est = garch_fit.ARCH{1};
                beta_est = garch_fit.GARCH{1};
                thetaM(2*i-1,1)=alpha_est;
                thetaM(2*i,1)=beta_est;
            end
                
            for t=1:T
                    Dt(:,:,t)=diag(cond_var(t,:));
            end
            
            prepared_data=std_returns;
                                              
    end
    
    % Devolver los datos preparados
    std_returns = prepared_data;
    
end 