%% Topoplotting mu_dks (sampled vs true)

handle = figure(30);
clf

cond = 1;

if cond == 1
    load('data/run_05bump_1000_1/histograms.mat');
else
    load('data/run_cond2_1/histograms.mat');
end

mu_dks = squeeze(hists.mu_dk(end, :, :));

annotation('textbox', [.0, .5, .5, .5], 'String', ['Sampled average Fan ' num2str(cond) ' bump topologies'], 'EdgeColor', 'None', 'FontSize', 22)

topo = reconstruct(squeeze(mu_dks(:, :)), coeff10, latent10, mean(data));
for k = 1:n
    xpos = (k - 1) / n -0.04;
    ypos = 0.65;

    axes('Units','normalized','Position',[xpos ypos .25 .25]);
    topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');
end

load('assoc_recog_classic_fan.mat')

annotation('textbox', [.0, .13, .5, .5], 'String', ['True Fan ' num2str(cond) ' bump topologies'], 'EdgeColor', 'None', 'FontSize', 22)

mags = magsb{cond};

topo = reconstruct(mags, coeff10, latent10, mean(data));
for k = 1:n
    xpos = (k - 1) / n -0.04;
    ypos = 0.3;

    axes('Units','normalized','Position',[xpos ypos .25 .25]);
    topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');

    if k == n
        colorbar('southoutside', 'Position', [0.035, 0.15, 0.9, 0.1], 'FontSize', 30)
    end
end

set(handle, 'Units', 'Inches');
pos = get(handle, 'Position');
set(handle, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

saveFile = fullfile('images', ['sampled_true_bumps_cond' num2str(cond)]);
% print(handle, saveFile, '-dpdf', '-r0');

%% Creating sampled and true stacked bar graph

n = 5;
cond = 2;

% Sampled proportions
if cond == 1
    handle = openfig('images/run_cond1/figure1.fig');
    
    set(handle, 'WindowStyle', 'normal', 'Position', [0, 0, 1400, 500])
    
    axObjs = handle.Children;
    subplot(1, 2, 1, axObjs(2))
    
    legend('off')
else
    handle = figure('Position', [0, 0, 1400, 600]);
    
    load('data/run_cond2_1/histograms.mat');
    t_is = squeeze(t_i_hist(end, :, :));
    
    subplot(1, 2, 1)

    barh(t_is, 'stacked')

    xlabel('t_i')
    ylabel('Participant')
end

title(['Stacked sampled Fan ' num2str(cond) ' duration proportions'])
set(gca, 'TitleHorizontalAlignment', 'left') 

set(gca, 'FontSize', 12)

annotation('textbox', [.08, .9, .1, .1], 'String', 'A', 'EdgeColor', 'None', 'FontSize', 22)

% True proportions
load('assoc_recog_overall_SNR5.mat');

lens = all_y - all_x + 1;
total_flat_lens = lens - n * 5;

bump_starts = all_pos;
bump_ends = bump_starts + 4;

flat_lens = [bump_starts zeros(length(lens), 1)] - 1;
flat_lens(:, 2:end-1) = bump_starts(:, 2:end) - bump_ends(:, 1:end - 1) - 1;
flat_lens(:, end) = lens - bump_ends(:, end);

flat_props = flat_lens ./ repmat(total_flat_lens, 1, n + 1);

nsubj = length(unique(all_subjects));
flat_props_subj = zeros(nsubj, n + 1);

for i = 1:nsubj
    flat_props_subj(i, :) = mean(flat_props(all_subjects == i & all_conds == cond, :), 1);
end

subplot(1, 2, 2)

barh(flat_props_subj, 'stacked')
legs = cell(n+1, 1);
for k = 1:n+1
    legs{k} = ['Flat ' num2str(k)];
end
set(gca, 'FontSize', 12)

title(['Stacked true Fan ' num2str(cond) ' duration proportions'])
set(gca, 'TitleHorizontalAlignment', 'left') 
xlabel('t_i')
ylabel('Participant')
xlim([0, 1])

legend(legs)

annotation('textbox', [.52, .9, .1, .1], 'String', 'B', 'EdgeColor', 'None', 'FontSize', 22)


set(handle, 'Units', 'Inches');
pos = get(handle, 'Position');
set(handle, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

saveFile = fullfile('images', ['sampled_true_proportions_cond' num2str(cond)]);
print(handle, saveFile, '-depsc', '-r0', '-painters');

%% Creating histograms for Fan 1 and Fan 2

n = 5;

burnIn = 1000;
thinLag = 20;
maxIter = 20000;

load('data/run_cond1/histograms.mat');
selection = 1:thinLag:maxIter-burnIn;
alpha1 = alpha_k_hist(selection, :);
beta1 = beta_k_hist(selection, :);

load('data/run_cond2_1/histograms.mat');
selection = burnIn+1:thinLag:maxIter;
alpha2 = alpha_k_hist(selection, :);
beta2 = beta_k_hist(selection, :);

handle = figure(4);
clf

for k = 1:n+1
    subplot(2, round((n+1)/2), k)
    hold on
    histogram(alpha1(:, k), 'FaceAlpha', 0.5, 'DisplayName', 'Fan 1')
    histogram(alpha2(:, k), 'FaceColor', '#D95319', 'FaceAlpha', 0.5, 'DisplayName', 'Fan 2')
    title(['\alpha_{' num2str(k) '}'])
    set(gca, 'FontSize', 14)
    
    if k == round((n+1)/2)
        legend
    end
    
    if mod(k-1, round((n+1)/2)) == 0
        ylabel('Counts')
    end
    
    if k > round((n+1)/2)
        xlabel('Parameter value')
    end
end

set(handle, 'Units', 'Inches');
pos = get(handle, 'Position');
set(handle, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

saveFile = fullfile('images', 'histograms_alpha_k');
print(handle, saveFile, '-dpdf', '-r0', '-painters');

handle = figure(5);
clf

for k = 1:n+1
    subplot(2, round((n+1)/2), k)
    hold on
    histogram(beta1(:, k), 'FaceAlpha', 0.5, 'DisplayName', 'Fan 1')
    histogram(beta2(:, k), 'FaceColor', '#D95319', 'FaceAlpha', 0.5, 'DisplayName', 'Fan 2')
    title(['\beta_{' num2str(k) '}'])
    set(gca, 'FontSize', 14)
    
    if k == round((n+1)/2)
        legend
    end
    
    if mod(k-1, round((n+1)/2)) == 0
        ylabel('Counts')
    end
    
    if k > round((n+1)/2)
        xlabel('Parameter value')
    end
end

set(handle, 'Units', 'Inches');
pos = get(handle, 'Position');
set(handle, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

saveFile = fullfile('images', 'histograms_beta_k');
print(handle, saveFile, '-dpdf', '-r0', '-painters');

%% Doing MAP inference for durations

n = 5;

% Cond 1 sampled

load('data/run_cond1/histograms.mat');

maxIter = 20000;
burnIn = 1000;
thinLag = 20;

displaySelection = 1:thinLag:maxIter-burnIn;

maps_modes = zeros(n + 1, 3);

for k = 1:n+1
    [counts, edges] = histcounts(alpha_k_hist(displaySelection, k));
    alpha_k_map = mean(edges(counts == max(counts)) + (edges(2) - edges(1)));
    
    [counts, edges] = histcounts(beta_k_hist(displaySelection, k));
    beta_k_map = mean(edges(counts == max(counts)) + (edges(2) - edges(1)));
    
    alpha_ik_mode = beta_k_map / (alpha_k_map + 1);
    
    maps_modes(k, :) = [alpha_k_map, beta_k_map, alpha_ik_mode];
end

alpha_iks = maps_modes(:, 3)';
duration_means1 = alpha_iks ./ sum(alpha_iks);

% Cond 2 sampled

load('data/run_cond2_1/histograms.mat');

maxIter = 20000;
burnIn = 1000;
thinLag = 20;

displaySelection = burnIn+1:thinLag:maxIter;

maps_modes = zeros(n + 1, 3);

for k = 1:n+1
    [counts, edges] = histcounts(alpha_k_hist(displaySelection, k));
    alpha_k_map = mean(edges(counts == max(counts)) + (edges(2) - edges(1)));
    
    [counts, edges] = histcounts(beta_k_hist(displaySelection, k));
    beta_k_map = mean(edges(counts == max(counts)) + (edges(2) - edges(1)));
    
    alpha_ik_mode = beta_k_map / (alpha_k_map + 1);
    
    maps_modes(k, :) = [alpha_k_map, beta_k_map, alpha_ik_mode];
end

alpha_iks = maps_modes(:, 3)';
duration_means2 = alpha_iks ./ sum(alpha_iks);

% True proportions
load('assoc_recog_overall_SNR5.mat');

lens = all_y - all_x + 1;
total_flat_lens = lens - n * 5;

bump_starts = all_pos;
bump_ends = bump_starts + 4;

flat_lens = [bump_starts zeros(length(lens), 1)] - 1;
flat_lens(:, 2:end-1) = bump_starts(:, 2:end) - bump_ends(:, 1:end - 1) - 1;
flat_lens(:, end) = lens - bump_ends(:, end);

flat_props = flat_lens ./ repmat(total_flat_lens, 1, n + 1);
flat_props1 = mean(flat_props(all_conds == 1, :), 1);
flat_props2 = mean(flat_props(all_conds == 2, :), 1);


h = figure('Position', [10, 10, 1000, 1000]);
clf

barh([duration_means2; flat_props2; nan(1, n+1); duration_means1; flat_props1], 'stacked')
title('Mean durations from across-subjects modes versus true durations')
xlabel('t_i')
yticklabels({'Sampled', 'True', '', 'Sampled', 'True'})
xlim([0, 1])

set(gca, 'FontSize', 17, 'TitleHorizontalAlignment', 'left')
annotation('textbox', [0, .64, .1, .1], 'String', 'Fan 1', 'EdgeColor', 'None', 'FontSize', 18, 'FontWeight', 'bold')
annotation('textbox', [0, .25, .1, .1], 'String', 'Fan 2', 'EdgeColor', 'None', 'FontSize', 18, 'FontWeight', 'bold')

legs = cell(n+1, 1);
for k = 1:n+1
    legs{k} = ['Flat ' num2str(k)];
end
legend(legs)

grid on

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

saveFile = fullfile('images', 'inferred_true_durations');
print(h, saveFile, '-dpdf', '-r0', '-painters');

%% Checking convergence

maxIter = 20000;
burnIn = 0;
thinLag = 1;

n = 5;

selection = burnIn+1:thinLag:maxIter;

runs = [1, 4];

for i = 1:length(runs)
    run = runs(i);
    load(['data/run_cond1_' num2str(run) '/histograms.mat'])
    
    
    % Flats
    figure(1)
    
    subplot(length(runs), 2, 2 * (i - 1) + 1);
    
    alpha_k_course = alpha_k_hist(selection, :);
    plot(selection, smoothdata(alpha_k_course, 'gaussian', 100), 'LineWidth', 2)
    
    title('\alpha_k')
    set(gca, 'FontSize', 16, 'TitleHorizontalAlignment', 'center')
    
    if 2 * (i - 1) + 1 <= (length(runs) - 1) * 2
        xticks([])
    else
        xlabel('Iteration')
    end
    
    if mod(2 * (i - 1), 2) == 0
        ylabel(['Run ' num2str(i)], 'FontWeight', 'bold')
    end
    
    subplot(length(runs), 2, 2 * i);

    beta_k_course = beta_k_hist(selection, :);
    plot(selection, smoothdata(beta_k_course, 'gaussian', 100), 'LineWidth', 2)
    
    title('\beta_k')
    set(gca, 'FontSize', 16, 'TitleHorizontalAlignment', 'center')
    
    if 2 * i <= (length(runs) - 1) * 2
        xticks([])
    else
        xlabel('Iteration')
    end
    
    if mod(2 * i - 1, 2) == 0
        ylabel(['Run ' num2str(i)], 'FontWeight', 'bold')
    end
    
    if 2 * i == 2 * length(runs)
        legs = cell(n+1, 1);
        for k = 1:n+1
            legs{k} = ['Flat ' num2str(k)];
        end
        legend(legs, 'Position', [0.93, 0.45, 0.05, 0.1])
    end

    
    % Bumps
    figure(2)

    for k = 1:n
        disp((i - 1) * length(runs) + k)
        subplot(length(runs), n, (i - 1) * n + k)
        
        mu_dk_course = mu_dk_hist(selection, :, k);
        plot(selection, smoothdata(mu_dk_course, 'gaussian', 100), 'LineWidth', 2)
        
        title(['\mu_{d' num2str(k) '}'])
        set(gca, 'FontSize', 14, 'TitleHorizontalAlignment', 'center')
        
        if mod((i - 1) * n + k - 1, n) == 0
            ylabel(['Run ' num2str(i)], 'FontWeight', 'bold')
        end
        
        if (i - 1) * n + k <= (length(runs) - 1) * n
            xticks([])
        else
            xlabel('Iteration')
            xticklabels({'0', '1e4', '2e4'})
        end
        
        if (i - 1) * n + k == n
            legs = cell(n+1, 1);
            for d = 1:10
                legs{d} = ['PC ' num2str(d)];
            end
            legend(legs, 'Position', [0.93, 0.47, 0.05, 0.1])
        end
    end
    
    h = figure(1);
    set(h, 'Units', 'Inches');
    pos = get(h, 'Position');
    set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
        'PaperSize', [pos(3), pos(4)])

    saveFile = fullfile('images', 'convergence_durations');
    print(h, saveFile, '-dpdf', '-r0', '-painters');
    
    h = figure(2);
    set(h, 'Units', 'Inches');
    pos = get(h, 'Position');
    set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
        'PaperSize', [pos(3), pos(4)])

    saveFile = fullfile('images', 'convergence_bumps');
    print(h, saveFile, '-dpdf', '-r0', '-painters');
end

%% Topoplotting mu_dks (fixed, known bumps)

handle = figure('Position', [0, 0, 1400, 550]);
clf

n = 5;

topo = reconstruct(squeeze(mu_dks(:, :)), coeff10, latent10, mean(data));
for k = 1:n
    xpos = (k - 1) / n - 0.02;
    ypos = 0.5;

    axes('Units','normalized','Position',[xpos ypos .25 .25]);
    topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');

    if k == n
        colorbar('southoutside', 'Position', [0.055, 0.2, 0.9, 0.1], 'FontSize', 30)
    end
end

set(handle, 'Units', 'Inches');
pos = get(handle, 'Position');
set(handle, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

saveFile = fullfile('images', 'average_topos_known_locs');
print(handle, saveFile, '-dpdf', '-r0');