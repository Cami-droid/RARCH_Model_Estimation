function [rotated_data, H_bar, Lambda, P] = rotate_data(outputs, model)
    T=outputs.T;
    returns=outputs.returns;
            
    H_bar = (1/T).*returns'*returns;
    % Descomposicion en valores propios
    [P, Lambda] = eig(H_bar);
    
    % Ensure that Lambda does not have elements close to zero to avoid numerical stability issues
    epsilon = 1e-10;
    Lambda_inv_sqrt = diag(1 ./ sqrt(max(diag(Lambda), epsilon)));

    % Calculate the inverse square root matrix of H_bar
    H_bar_inv_sqrt = P * Lambda_inv_sqrt * P';

    
    switch model
        case {'RBEKK'}

            % Rotation for RBEKK y GOGARCH
            e_t = returns*H_bar_inv_sqrt;

             
        case {'OGARCH'}

            % Specific rotation for OGARCH;
            e_t = returns*H_bar_inv_sqrt * P';

        case {'GOGARCH'}
            
            % Specific rotation for RBEKK and GOGARCH
            e_t = returns*H_bar_inv_sqrt;
            
        case {'RDCC'}
            PI_bar=H_bar; %%  returns tha income from prepare_data are actually standarized for RDCC
            %%% i'm just giving it the names that are on the paper to avoid confusions
            
                        
            % specific rotation for RDCC
            PI_bar_inv_sqrt = P * sqrt(inv(Lambda)) * P';        
            e_t = returns*PI_bar_inv_sqrt;
            
       
    end
    
    rotated_data = e_t;
end
