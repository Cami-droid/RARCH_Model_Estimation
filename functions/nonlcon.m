% Non linear restrictions function for Common Persistence specifiction

function [c, ceq] = nonlcon(params)

    alpha   =   params(1);
    beta    =   params(2);
    lambda  =   params(3);
    
    % Restricciones de desigualdad (c <= 0)
    c = [max(alpha) - lambda; % \lambda >= max(\alpha_{ii})
         -lambda;             % \lambda > 0, se maneja con lb
         lambda - 1];         % \lambda < 1, se maneja con ub y b
    
    %% There isn't any equality restrictions (ceq = []) %%
    ceq = [];
end