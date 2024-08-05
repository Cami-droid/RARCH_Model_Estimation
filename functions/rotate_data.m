function [rotated_data, H_bar, Lambda, P ,PI_bar, Lambda_C,P_C] = rotate_data(outputs, model)
    T=outputs.T;
    rt=outputs.returns';
    Et=outputs.std_returns';
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% H_bar Eigenvalues Decomposition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_bar = (1/T).*rt*rt';
    [P, Lambda] = eig(H_bar);
    %P=[0.394,-0.919;0.919,0.394];  % these are the Table 1 values in the paper for P
    %Lambda=diag([1.820,9.882]);    % these are the Table 1 values in the paper for Lambda
    epsilon = 1e-10;
    % Lambda=diag(max(diag(Lambda),epsilon));
    % Ensure that Lambda does not have elements close to zero to avoid numerical stability issues
    Lambda_inv_sqrt = diag(1 ./ sqrt(diag(Lambda)));
    % Calculate the inverse square root matrix of H_bar
    H_bar_inv_sqrt = P * Lambda_inv_sqrt * P';
    PI_bar=NaN(size(H_bar));
    Lambda_C=NaN(size(Lambda));
    P_C=NaN(size(P));

    switch model
        case {'RBEKK'}

            % Rotation for RBEKK y GOGARCH
            e_t = H_bar_inv_sqrt*rt;
            
             
        case {'OGARCH'}
            PI_bar=[];
            Lambda_C=[];
            P_C=[];
            % Specific rotation for OGARCH;
            e_t = Lambda_inv_sqrt*P'*rt;

        case {'GOGARCH'}
            
            % Specific rotation for RBEKK and GOGARCH
            e_t = H_bar_inv_sqrt*rt;
            
        case {'RDCC'}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PI_bar Eigenvalues Decomposition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Et is a more accurate name, they are standarized returns;
        PI_bar=(1/T).*Et*Et';

        [P_C,Lambda_C]=eig(PI_bar); %%  returns tha income from prepare_data are actually standarized for RDCC
        %P_C=[0.706,-0.708;0.708,0.706]; % these are the Table 1 values in the paper for P_C
        %Lambda_C=diag([0.512,1.452]);   % these are the Table 1 values in the paper for Lambda_C
        %Lambda_C=diag(max(diag(Lambda),epsilon));

        % Ensure that Lambda does not have elements close to zero to avoid numerical stability issues

        Lambda_C_inv_sqrt = diag(1 ./ sqrt(diag(Lambda_C)));
        PI_bar_inv_sqrt = P_C * Lambda_C_inv_sqrt * P_C'; 

            
            e_t = PI_bar_inv_sqrt*Et; % e_t stands for rotated returns
                
    end
    
    rotated_data = e_t';
end
