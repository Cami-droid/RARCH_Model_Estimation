function run_rarch_models(table_number)
    % Run RARCH models based on the specified table number
    % 
    % Inputs:
    %   table_number - The number of the table to be processed (2, 4, or 5)
    
    % Define the relative path to the results directory
    results_dir = fullfile(pwd, 'results');
    
    % Open log file
    logFilePath = fullfile(results_dir, 'run_rarch_models_log.txt');
    diary(logFilePath);
       
    % Validate the input table number
    if ~ismember(table_number, [2, 4, 5])
        error('Invalid table number. Please enter 2, 4, or 5.');
    end
    
    % Define the name of the data import file based on the selected table
    data_import_file = sprintf('data_import_Table%d', table_number);

    % Check if the data import file exists
    if exist([data_import_file, '.m'], 'file') == 2
        % Execute the data import file
        run(data_import_file);
        
        % Execute the main script based on the table number
        if table_number == 5
            main_5;
        else
           main;
        end

        % Save all variables to a .mat file
        results_dir = fullfile(pwd, 'results');
        save(fullfile(results_dir, sprintf('results_table%d.mat', table_number)));
    else
        error('The file %s does not exist.', data_import_file);
    end
    
    % Close the log file
    diary off;
end
