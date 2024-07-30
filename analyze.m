files = dir('**/*.m');
for k = 1:length(files)
    results = checkcode(fullfile(files(k).folder, files(k).name));
    if ~isempty(results)
        fprintf('Issues in %s:\n', files(k).name);
        disp(results);
    end
end
