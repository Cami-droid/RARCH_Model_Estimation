function generate_table(data)
    % This function generates a table from the provided structure.
    % data - structure with fields: model, specification, params, cond_cov

    % Extract fields from the structure
    model = data.model;
    specification = data.specification;
    params = data.params;
%     cond_cov = data.cond_cov;
    
    % Display the table headers
    fprintf('Model: %s\n', model);
    fprintf('Specification: %s\n', specification);
    
    % Display the parameters
    fprintf('Parameters:\n');
    disp(params);
    
    %% Display the conditional covariance
    %fprintf('Conditional Covariance:\n');
%     disp(cond_cov); %%
    
    % Example of creating a table using MATLAB's table function
    % Assuming params and cond_cov are numeric arrays
%     param_table = array2table(params, 'VariableNames', {'Param1', 'Param2', 'Param3'}); % Modify variable names as needed
%     cond_cov_table = array2table(cond_cov, 'VariableNames', {'Cov1', 'Cov2', 'Cov3'}); % Modify variable names as needed
    
    % Combine tables (if desired)
%     combined_table = [param_table, cond_cov_table];
    
    % Display the combined table
%     disp(param_table);
    
    % Save the table to a file (optional)
%     writetable(param_table, 'results_table.csv'); % Save as CSV
end
