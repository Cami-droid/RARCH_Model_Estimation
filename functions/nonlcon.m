% Non linear restrictions function for Common Persistence specification
function [c, ceq] = nonlcon(specification,,outputs,thetaD)
    d=outputs.d

    switch specification
        case 'Scalar'
        
        alpha = thetaD(1) ;
        beta = thetaD(2);

        c = [-alpha; %alpha>=0 (has to be strictly positive)
            - beta]   %beta>=0                       
                    
    
        case 'Diagonal'

        alpha = thetaD(1:d);
        beta = thetaD(d+1:2*d);

        c1=-eye(2*d)  %alphas and betas >=0
        c2=eye(d),eye(d)

        case 'CP'

        alpha_vec = thetaD(1:d) ;
        lambda_cp = thetaD(d+1) ;
        beta_vec =(lambda_cp*ones(1,d) - alpha_vec); %% thetaD(d+1) is lambda and the thetaD(1:d) are alphas, all term root squared
     
        c = [max(alpha) - lambda;                           % \lambda >= max(\alpha_{ii})
                            -lambda;             % \lambda > 0, as lb
                            lambda - 1;          % \lambda < 1, as ub y b
                -alpha                 ];             % alpha > 0

    end
    
    % Inequality restrictions (c <= 0)
    c = [max(alpha_vec) - lambda; % \lambda >= max(\alpha_{ii})
         -lambda;             % \lambda > 0, as lb
         lambda - 1;          % \lambda < 1, as ub y b
         -alpha];             % alpha > 0
    
    %% There isn't any equality restrictions (ceq = []) %%
    ceq = [];
end
