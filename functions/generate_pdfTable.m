% Create a figure and plot the table
f = figure('Position', [100, 100, 1200, 400]); % Adjust the figure size as needed

% Create a UITable within the figure
t = uitable('Data', thetaD_table_matlab.Variables, ...
            'ColumnName', thetaD_table_matlab.Properties.VariableNames, ...
            'RowName', thetaD_table_matlab.Properties.RowNames, ...
            'Units', 'Normalized', ...
            'Position', [0, 0, 1, 1]);

% Adjust the table size to fit the data
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

% Adjust the figure size to match the table size
f.Position(3) = t.Position(3);
f.Position(4) = t.Position(4);

% Set the background color of the figure to white (to ensure visibility)
f.Color = 'r';

% Define the directory and file name for the PDF
results_dir = 'D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\results';
pdf_file = fullfile(results_dir, 'thetaD_table.pdf');

% Save the figure as a PDF
print(f, pdf_file, '-dpdf', '-bestfit');

% Optionally display the figure to check the result
uiwait(f); % Uncomment to display the figure

% Close the figure
close(f);
