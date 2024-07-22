function Gt = calc_all_Gts(model, specification, outputs, thetaD, initial_Gt)
% the initial subindex (which in the theory is 0, is 1 in this function)
T=outputs.T;
d=outputs.d;
rotated_returns=outputs.rotated_returns;
thetaS=outputs.H_bar;

Gt(:,:,1)=initial_Gt;

for t_count= 2:T+1 % Adjusting loop to account for t starting from 0
        
        Gt(:, :, t_count) = calcGt(model,specification,outputs, thetaD, Gt(:, :, t_count-1)); % t_count-1=t for returns to match theory
        
        % Adding a small regularization term to Gt to avoid singularity
        reg_term = 1e-6 * eye(d);
        Gt_reg = Gt(:, :, t_count) + reg_term;

end