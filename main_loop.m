i = 1;
settings(i).path = '~/Documents/Uni/FYRP/';
settings(i).save_dir = 'data/run_1bump';
settings(i).dataset = "";
settings(i).visualize_flag = false;
settings(i).parameters.n = 1;
settings(i).parameters.N = 20;

i = i + 1;
settings(i).path = '~/Documents/Uni/FYRP/';
settings(i).save_dir = 'data/run_2bump';
settings(i).dataset = "";
settings(i).visualize_flag = false;
settings(i).parameters.n = 2;
settings(i).parameters.N = 20;

p = gcp('nocreate');
if isempty(p)
    parpool(7)
end

parfor i = 1:2
    s = settings(i);
    main(s.path, s.save_dir, s.dataset, s.visualize_flag, s.parameters);
end