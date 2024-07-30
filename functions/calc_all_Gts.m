function Gt = calc_all_Gts(model, specification, outputs, thetaD, initial_Gt)
% the initial subindex (which in the theory is 0, is 1 in this function)
T=outputs.T;
d=outputs.d;
rotated_returns=outputs.rotated_returns;
thetaS=outputs.H_bar;

Gt(:,:,1) = calcGt(model,specification,outputs, thetaD, initial_Gt,1);

for t= 2:T % Adjusting loop to account for t starting from 1
        
        Gt(:, :, t) = calcGt(model,specification,outputs, thetaD, Gt(:, :, t-1),t); 
        
        % Adding a small regularization term to Gt to avoid singularity
        reg_term = 1e-6 * eye(d);
        Gt(:,:,t) = Gt(:,:,t) + reg_term;
end

