
function [parameters,ll,Ht,VCV,scores]=gogarch_spec(data,p,q,gjrType,type,startingVals,options,specification)
%(data,p,q,gjrType,type,startingVals,options
% the same gogarch function in mfe toolkit with the specification as additional input
% Specification could be 'Scalar','Diagonal','CP' which is common persistence. 'Diagonal' specification returns 
% the same output as the original gogarch function.
% type:'OGARCH'or 'GOGARCH'
d=size(data,2);
idxS= (d * (d + 1)) / 2;
idxPhis = (d * (d - 1)) / 2;

param_Phis=NaN(idxPhis,1);
param_thetaS=NaN(idxS,1);
param_delta=NaN(1,1);

        

method='2-stage';

    switch specification

        case 'Scalar'

        param_alphas=NaN(1,1); %momentaneo
        param_betas=NaN(1,1);
                
                
        switch type
            case 'GOGARCH'
                parameters=[param_thetaS;param_Phis;param_alphas;param_betas;param_delta];
            case 'OGARCH'
                parameters=[param_thetaS;param_alphas;param_betas];
        end

        ll=NaN;
        Ht=NaN;
        VCV=NaN;
        scores=NaN;
  

        case 'Diagonal'
        
        [parameters,ll,Ht,VCV,scores] = gogarch(data,p,q,gjrType,type,startingVals,options);

        switch type
            case 'GOGARCH'
                parameters=[parameters;param_delta];
            case 'OGARCH'
                parameters;
        end

        case 'CP'

            
        param_alphas=NaN(d,1); %momentaneo
        param_lambda_cp=NaN(1,1);
        
        
        switch type
            case 'GOGARCH'
            
            parameters=[param_thetaS;param_Phis;param_alphas;param_lambda_cp;param_delta];

            case 'OGARCH'
            
            parameters=[param_thetaS;param_alphas;param_lambda_cp];
        
        end

        ll=NaN;
        Ht=NaN;
        VCV=NaN;
        scores=NaN;
       

        otherwise

        error('Specification not supported');

    end

    parameters=parameters';
    
end




                