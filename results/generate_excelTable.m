% Define the filename and path for the Excel file
results_dir = 'D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\results';
excel_file = fullfile(results_dir, 'thetaD_table.xlsx');

% Write the MATLAB table to an Excel file
writetable(thetaD_table_matlab, excel_file, 'WriteRowNames', true);

% Display a message indicating success
disp(['Table has been successfully exported to ', excel_file]);
