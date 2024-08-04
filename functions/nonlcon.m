% Non linear restrictions function for Common Persistence specification
function [c, ceq] = nonlcon(params)

    alpha   =   params(1);
    beta    =   params(2);
    lambda  =   params(3);
    
    % Inequality restrictions (c <= 0)
    c = [max(alpha) - lambda; % \lambda >= max(\alpha_{ii})
         -lambda;             % \lambda > 0, as lb
         lambda - 1;          % \lambda < 1, as ub y b
         -alpha];             % alpha > 0
    
    %% There isn't any equality restrictions (ceq = []) %%
    ceq = [];
end
