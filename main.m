% Gibbs sampler for discovering cognitive stages in EEG data

% There are many variables that need to be sampled in turn:
% - We start the entire procedure by fixing every variable to some non-
% specific value.
% - Then, repeat until convergence:
%   - Sample every variable within the flat box. For every flat, we
%   sample alpha_k and beta_k. Then, for every subject, sample alpha_ik.
%   - Consequently, sample the remaining flat variable for every subject, 
%   namely t_i.
%   - Finally, sample mu_dk and tau^2_dk for every bump and every principal
%   component. Additionally, sample m_idk for every subject.

%% Loading data

path = '~/Documents/Uni/FYRP/';
cd(path);

addpath([path 'data'], ... 
    genpath([path 'synth_data']), ...
    genpath([path 'hierarchical_gibbs']), ...
    [path 'eeglab2021.0']);

% A number will automatically be appended to this path to be able to
% discern between different runs with the same savePath (e.g. 'data/run'
% will become 'data/run_1', after that 'data/run_2', etc.)
savePath = [path 'data/run_fixed_known'];

load('varForBumps135_100.mat');
load('assoc_recog_overall_SNR5.mat');
load('chanlocs.mat');

x = all_x;
y = all_y;
subjects = all_subjects;
conds = all_conds;
normedscore10 = all_signal;

trial_lens = y - x + 1;

%% Setting up EEGlab

eeglab;
close % the EEGlab figure

%% Setting up variables

max_iter = 20000;

% Which condition to sample for.
cond = 1;

% Number of bumps, so there are n+1 flats.
n = 5;
% A normalising variable from Anderson et al. (2016)
V = 5;

N = length(unique(subjects));
% N = 2;
D = size(normedscore10, 2);

% Bump topology - 1 by 5 matrix
P = sin(pi/10:pi/5:pi-pi/10);

% Note: The appended 's' in a variable name indicates multiple instances of 
% that variable

%%% Hyperparameters, will not be sampled
% Set so the hyperpriors are vague

% Flat box
a_k = 1;
b_k = 0.5;
q_k = 1;
c_k = 2;
d_k = 0.1; % Rate parameter

% Bump box
f_dk = 0;
g2_dk = 1000;
l_dk = 1;
o_dk = 0.01;

%%% Parameters, will be sampled

% Flat box
alpha_ks = ones(n+1, 1);
beta_ks = ones(n+1, 1);
alpha_iks = zeros(N, n+1); alpha_iks(:) = 2;

% Last variable on flat branch
% Initially assume equally distributed flats
t_is = zeros(N, n+1); t_is(:) = 1 / (n+1);

% Bump box
mu_dks = zeros(D, n); % Principal component box
tau2_dks = ones(D, n); % Principal component box
m_idks = zeros(N, D, n);

%% Setting up parameter histogram variables

alpha_k_hist = zeros(max_iter, n+1);
beta_k_hist = zeros(max_iter, n+1);
alpha_ik_hist = zeros(max_iter, N, n+1);
t_i_hist = zeros(max_iter, N, n+1);
m_idk_hist = zeros(max_iter, N, D, n);
mu_dk_hist = zeros(max_iter, D, n);
tau2_dk_hist = zeros(max_iter, D, n);

%% Sampling

tic

for iter = 1:max_iter
    disp(iter)
    
    % Record parameter values into histograms variables
    alpha_k_hist(iter, :) = alpha_ks;
    beta_k_hist(iter, :) = beta_ks;
    alpha_ik_hist(iter, :, :) = alpha_iks;
    t_i_hist(iter, :, :) = t_is;
    m_idk_hist(iter, :, :, :) = m_idks;
    mu_dk_hist(iter, :, :) = mu_dks;
    tau2_dk_hist(iter, :, :) = tau2_dks;
    
    % Visualize results
    visualize(max_iter, iter, ... 
        normedscore10, x, y, subjects, ... 
        coeff10, latent10, data, chanlocs, ...
        alpha_k_hist, beta_k_hist, alpha_ik_hist, t_i_hist, ... 
        mu_dk_hist, tau2_dk_hist, m_idk_hist);
    
    disp('Sampling flat box...')
    for k = 1:n+1
        % For every flat, we sample alpha_k and beta_k.
        alpha_ks(k) = alpha_kMH(alpha_ks(k), a_k, b_k, q_k, alpha_iks(:, k), beta_ks(k), 1);
        beta_ks(k) = gamrnd(N * a_k + c_k, 1/(sum(1 ./ alpha_iks(:, k)) + d_k));
        
        for i = 1:N
             % Then, for every subject, sample alpha_ik.
             alpha_iks(i, k) = alpha_ikMH(alpha_iks(i, setdiff(1:n+1, k)), alpha_iks(i, k), t_is(i, k), alpha_ks(k), beta_ks(k), 1);
        end
    end
    
    % Consequently, sample the remaining flat variable for every subject, 
    % namely t_i.
    disp('Sampling durations...')
    for i = 1:N
        t_is(i, :) = tiMH(...
            normedscore10, t_is(i, :), alpha_iks(i, :), V, ... 
            x(subjects == i & conds == cond), ... 
            y(subjects == i & conds == cond), n, 1);
    end
    
    disp('Sampling bump box...')
    for d = 1:D
        for k = 1:n
            % Finally, sample mu_dk and tau^2_dk for every bump and every
            % principal component.
            mu_dks(d, k) = normrnd( ... 
                (g2_dk * sum(m_idks(:, d, k)) + tau2_dks(d, k) * f_dk) / (N * g2_dk + tau2_dks(d, k)), ...
                (tau2_dks(d, k) * g2_dk) / (N * g2_dk + tau2_dks(d, k)) );
            tau2_dks(d, k) = invgamrnd(N / 2 + l_dk, ...
                (sum((m_idks(:, d, k) - mu_dks(d, k)).^2) + 2 * o_dk) / 2);
            
            % Additionally, sample m_idk for every subject.
            for i = 1:N
                % This requires some work. We need to determine where every
                % bump is located given t_is - this is different for every
                % trial, as t_is only contains proportions.
                
                % First, acquire all data for this subject.
                subject_x = x(subjects == i & conds == cond);
                subject_lens = trial_lens(subjects == i & conds == cond);
                
                bump_indices = getBumpIndices(subject_x, subject_lens, ...
                    t_is(i, :), k, n);
                
                % C: Number of trials
                C = length(subject_x);
                
                % Due to rounding errors, it is possible that the bump
                % indices extend just beyond the data array boundaries;
                % correct for that by padding.
                bump_indices = bump_indices(bump_indices < length(normedscore10));
                pad_size = length(P) * C - length(bump_indices);

                % Get those data!
                bumps = normedscore10(bump_indices, d);
                bumps = padarray(bumps, [pad_size 0], 'post');
                
                bump_corr = repmat(P, 1, C) * bumps;
                
                % And finally draw the sample.
                denominator = 2 * tau2_dks(d, k) * C * sum(P .^ 2) + V;
                m_idks(i, d, k) = normrnd( ... 
                    (2 * tau2_dks(d, k) * bump_corr + V * mu_dks(d, k)) / denominator, ...
                    (V * tau2_dks(d, k)) / denominator);
            end
        end
    end
end

toc

%% Save histogram variables to file

% Determine which number to append to the save directory
dirs = regexp(savePath, filesep, 'split');
parentDir = fullfile(dirs{1:end-1});
saveDir = dirs{end};

pattern = ['^' saveDir '_?([0-9]+)?$'];

parentDirItems = dir(parentDir);

maxSaveDirNum = 0;
for i = 1:length(parentDirItems)
    token = regexp(parentDirItems(i).name, pattern, 'tokens');
    if ~isempty(token)
        saveDirNumber = str2double(token{1, 1});
        if ~isnan(saveDirNumber)
            maxSaveDirNum = max(saveDirNumber, maxSaveDirNum);
        end
    end
end

saveDir = [saveDir '_' num2str(maxSaveDirNum + 1)];

savePath = fullfile(parentDir, saveDir);
mkdir(savePath);

save(fullfile(savePath, 'histograms'), ... 
    'alpha_k_hist', ...
    'beta_k_hist', ...
    'alpha_ik_hist', ...
    't_i_hist', ...
    'm_idk_hist', ...
    'mu_dk_hist', ...
    'tau2_dk_hist');
