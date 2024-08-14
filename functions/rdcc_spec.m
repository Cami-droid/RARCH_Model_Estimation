function [parameters, ll ,Ht, VCV, scores, diagnostics]=rdcc_spec(data,dataAsym,m,l,n,p,o,q,gjrType,method,composite,startingVals,options,specification)
d=size(data,2);
% the same dcc function in mfe toolkit with the specification as additional input 
% and with a rotation of returns  after its standarization

% Specification could be 'Scalar','Diagonal','CP' which is common persistence. 'Diagonal' specification returns 
% the same output as the original gogarch function.
method='2-stage';
idxM = 3 * d; 
idxS = (d * (d - 1)) / 2; %correlations, diagonal is excluded

    switch specification
        case 'Scalar'

        [parameters, ll ,Ht, VCV, scores, diagnostics]=dcc(data,dataAsym,m,l,n,p,o,q,gjrType,method,composite,startingVals,options);
    

        case 'Diagonal'

        
        idxD = 2*d;

        parameters=NaN(1,idxM+idxS+idxD); %momentaneo
        ll=NaN;
        Ht=NaN;
        VCV=NaN;
        scores=NaN;
        diagnostics=NaN;
        
        
        case 'CP'
        
        idxD = d+1;

        parameters=NaN(1,idxM+idxS+idxD); %momentaneo 
        ll=NaN;
        Ht=NaN;
        VCV=NaN;
        scores=NaN;
        diagnostics=NaN;
                 
        otherwise

        error('Specification not supported');

    end

    
        
end