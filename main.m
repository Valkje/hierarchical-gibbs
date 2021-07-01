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

addpath([path 'data'], genpath([path 'synth_data']), ...
    genpath([path 'hierarchical_gibbs']), ...
    [path 'eeglab2021.0']);

load('varForBumps135_100.mat');
load('assoc_recog_overall_SNR5_lockedPos.mat');
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

% Which condition to sample for.
cond = 2;
pos = all_pos(conds == cond, :);

% Number of bumps, so there are n+1 flats.
n = 5;
% A normalising variable from Anderson et al. (2016)
V = 5;

% N = length(unique(subjects));
N = 2;
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

%% Calculating true t_is

true_n = size(all_pos, 2);
num_trials = length(x);

true_bump_ends = zeros(num_trials, true_n+1);
true_bump_ends(:, 2:end) = all_pos + 4;
true_bump_starts = zeros(num_trials, true_n+1) + trial_lens + 1;
true_bump_starts(:, 1:true_n) = all_pos;
true_t_is = true_bump_starts - true_bump_ends - 1; % In number of samples

total_flat_lens = trial_lens - true_n * 5;
true_t_is = true_t_is ./ repmat(total_flat_lens, 1, true_n + 1);
true_t_is = mean(true_t_is, 1); % In average proportions

%% Sampling

max_iter = 20000;

for iter = 1:max_iter
    disp(iter)
    visualize(max_iter, iter, ... 
        normedscore10, x, y, subjects, ... 
        coeff10, latent10, data, chanlocs, ...
        alpha_ks, beta_ks, alpha_iks, t_is, mu_dks, tau2_dks, m_idks);
    
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
%         t_is(i, :) = true_t_is;
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

%                 bump_starts = repmat(pos(1, k), C, 1);
%                 
%                 % To index into the entire data array, we need to offset
%                 % the bump locations by every trial's onset.
%                 bump_starts = subject_x + bump_starts - 1; % C by 1
%                 
%                 % Determine bump indices, reshape into column vector.
%                 bump_indices = repmat(bump_starts, 1, 5) + (0:4);
%                 bump_indices = reshape(bump_indices', [], 1);
                
                % Get those data!
                bumps = normedscore10(bump_indices, d);
                
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