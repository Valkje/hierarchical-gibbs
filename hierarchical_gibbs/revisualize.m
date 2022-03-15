% Script used for most of the images that are in the extended report. For
% many subjects and bumps, it is awfully slow (i.e. it can take a minute or
% two before all figures have loaded).

path = '~/Documents/Uni/FYRP/';
hists_file = fullfile(path, ...
    'data/run_N20_n10_noSignHack_1/histograms.mat');

cd(path);

addpath(genpath([path 'data']), ... 
    genpath([path 'synth_data']), ...
    genpath([path 'hierarchical_gibbs']), ...
    [path 'eeglab2021.0']);

load varForBumps135_100.mat coeff10 latent10 data
load chanlocs.mat chanlocs

if contains(hists_file, 'synthetic')
    load assoc_recog_overall_SNR5.mat all_x all_y all_subjects all_conds ... 
        all_signal all_pos

    x = all_x;
    y = all_y;
    subjects = all_subjects;
    conds = all_conds;
    normedscore10 = all_signal;
else
    load varForBumps135_100.mat x y subjects conds normedscore10
end

% eeglab;
% close

%% Loading hists and initialising visualisations

load(hists_file);

params.max_iter = size(hists.xi_k, 1);
params.iter = params.max_iter;
params.start_iter = 1;
params.normedscore10 = normedscore10;
params.x = x;
params.y = y;
params.subjects = subjects;
params.coeff10 = coeff10;
params.latent10 = latent10;
params.data = data;
params.chanlocs = chanlocs;
params.conds = conds;
params.cond = 1;
params.relative = false;
params.firstVis = true;

params.burnIn = 10;
params.renderLag = 1;
params.thinLag = 1;

visualize(params, hists);

%% Adding colorbar to topologies, save file

save_path = '~/Documents/Uni/FYRP/images/extended_report';

fig = figure(10);
% set(fig, 'WindowStyle', 'normal', 'Units', 'pixels', 'Position', [0 0 1200 800])

axes('Units','normalized','Position',[0.9 0.03 0.1 0.95])
axis off
% rectangle(ax, 'Position', [0.93, 0.05, 0.135, 0.91], 'Curvature', 0.2, 'FaceColor', [1 1 1])
rectangle('Position', [0, 0, 1, 1], 'FaceColor', [1 1 1, 0.5])

colorbar('off')
colorbar('eastoutside', 'Position', [0.93, 0.05, 0.025, 0.9], 'FontSize', 20, 'Limits', [-12 12], 'AxisLocation', 'in')
caxis([-12 12])

path_parts = split(hists_file, filesep);
save_file = fullfile(save_path, path_parts{end-1});

savePDF(save_file)
saveas(gcf, [save_file '.png'])

%% Topoplotting mu

figure(24)
clf

% Sampled
mu_dks = squeeze(hists.mu_dk(params.iter, :, :));
n = size(mu_dks, 2);
ax_size = min(1 / n, 0.4);

annotation('textbox', [.0, .5, .5, .5], 'String', ['Sampled bump topologies'], 'EdgeColor', 'None', 'FontSize', 22)

topo = reconstruct(squeeze(mu_dks(:, :)), coeff10, latent10, mean(data));
for k = 1:n
    xpos = (k - 1) / n;
    ypos = 0.7;
    
    axes('Units','normalized','Position',[xpos ypos ax_size ax_size]);
    topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');
end

% True
load('assoc_recog_classic_fan.mat')

ax_size = 1/5;

annotation('textbox', [.0, .13, .5, .5], 'String', ['True bump topologies'], 'EdgeColor', 'None', 'FontSize', 22)

mags = magsb{1}; % cond == 1

topo = reconstruct(mags, coeff10, latent10, mean(data));
for k = 1:5
    xpos = (k - 1) / 5;
    ypos = 0.3;
    
    axes('Units','normalized','Position',[xpos ypos ax_size ax_size]);
    topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');

    if k == 5
        colorbar('southoutside', 'Position', [0.05, 0.15, 0.9, 0.1], 'FontSize', 30)
    end
end

% Save figure as png
saveas(gcf, [save_file '_mu.png'])