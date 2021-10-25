fileID = fopen('settings.json', 'r');
contents = fscanf(fileID, '%s');
fclose(fileID);

contents = jsondecode(contents);
settings = contents.settings;

p = gcp('nocreate');
if isempty(p)
    parpool(7)
end

parfor i = 1:length(settings)
    s = settings(i);
    disp(s.parameters.n)
%     main(s.path, s.save_dir, s.dataset, s.visualize_flag, s.parameters);
end