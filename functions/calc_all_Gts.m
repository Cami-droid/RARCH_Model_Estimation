function Gt = calc_all_Gts(thetaS, thetaD, rotated_returns, initial_Gt, model, specification)
% the initial subindex (which in the theory is 0, is 1 in this function)
T=size(rotated_returns,1)
d=size(rotated_returns,2)
Gt(:,:,1)=initial_Gt;

for t = 2:T+1 % Adjusting loop to account for t starting from 0
        Gt(:, :, t) = calcGt(thetaS, thetaD, rotated_returns(t-1, :), Gt(:, :, t-1), model, specification); % t-1 for returns to match theory
        
        % Adding a small regularization term to Gt to avoid singularity
        reg_term = 1e-6 * eye(d);
        Gt_reg = Gt(:, :, t) + reg_term;

end