function [] = main(path, save_dir, dataset, visualize_flag, parameters)
% Gibbs sampler for discovering cognitive stages in EEG data.
% Arguments:
% - path: Path to the folder that contains this file. That folder should
%   contain the subfolders 'data', 'synth_data', 'hierarchical_gibbs', and
%   'eeglab2021.0'.
% - save_dir: Directory in which the values of all parameters, for every
%   time step, will be saved. On saving, this directory will automatically
%   be appended with a version number to prevent overwriting old data.
% - dataset: The synthetic dataset to use for estimation. If the estimation
%   should be done with human data instead, specify this as an empty
%   string.
% - visualize_flag: Whether to create intermediate visualizations while
%   sampling. Especially for many bumps and participants, turning this off
%   can be beneficial.
% - parameters: A struct that specifies some parameters for the algorithm.
%   Currently, the only supported fields are 'n' (the number of bumps) and
%   'N' (the number of participants).

% There are many variables that need to be sampled in turn:
% - We start the entire procedure by fixing every variable to some non-
%   specific value. (Starting values should, in principle, not matter.)
% - Then, repeat until convergence:
%   - Sample every variable within the left branch. For every bump, we sample
%     xi_k, upsilon2_k, alpha_k, and beta_k. Then, for every subject, sample
%     nu_ik and sigma2_ik. Finally, for every trial, sample p_itk. 
%   - As for the right branch, sample mu_dk and tau^2_dk for every
%     bump and every principal component. Additionally, sample m_idk for
%     every subject.

arguments
    path char = '~/Documents/Uni/FYRP/'
    save_dir char = 'data/run_debug'
    dataset string = ""
    visualize_flag logical = false
    parameters struct = struct
end

%% Checking validity of and setting defaults for PARAMETERS

if ~isfield(parameters, 'n')
    parameters.n = 5; % Number of bumps
end

if ~isfield(parameters, 'N')
    parameters.N = 20; % Number of subjects
end

%% Loading data

cd(path);

addpath(genpath([path 'data']), ... 
    genpath([path 'synth_data']), ...
    genpath([path 'hierarchical_gibbs']), ...
    [path 'eeglab2021.0']);

% A number will automatically be appended to this path to be able to
% discern between different runs with the same savePath (e.g. 'data/run'
% will become 'data/run_1', after that 'data/run_2', etc.)
savePath = [path save_dir];

load varForBumps135_100.mat coeff10 latent10 conds data normedscore10 subjects x y
load chanlocs.mat chanlocs

% If DATASET is filled, assume that we want to run the algorithm with a
% synthetic dataset.
if ~(dataset == "")
    load(dataset, 'all_x', 'all_y', 'all_subjects', 'all_conds', ... 
        'all_signal')
    
    x = all_x;
    y = all_y;
    subjects = all_subjects;
    conds = all_conds;
    normedscore10 = all_signal;
end

trial_lens = y - x + 1;

% Whether to continue with a previous run. If undesired, please comment
% this line.
continue_path = [path 'data/run_N20_n10_incremental_mode_1'];

%% Setting up EEGlab

if visualize_flag
    eeglab;
    close % the EEGlab figure
end

%% Setting up variables

max_iter = 500;

% Used to initialize figures
firstVis = true;

% Which condition to sample for.
cond = 1;

% Number of bumps, so there are n+1 flats.
n = parameters.n;
% A normalising variable from Anderson et al. (2016)
V = 5;

% N = length(unique(subjects));
N = parameters.N;
D = size(normedscore10, 2);

% Bump topology - 1 by 5 matrix
P = sin(pi/10:pi/5:pi-pi/10);

% Note: The appended 's' in a variable name indicates multiple instances of 
% that variable

%%% Hyperparameters, will not be sampled
% Set so the hyperpriors are vague

% Left branch
% xi_k prior (normal)
n_ks = ((1:n) * (mean(trial_lens) / (n + 1)))'; % Even spread of bumps as prior
% n_ks =  all_pos(1, :)';
s2_k = 100;
% upsilon2_k prior (inverse gamma)
a_k = 1;
b_k = 0.01;
% alpha_k prior (non-closed form)
c_k = 50; 
d_k = 0.5;
e_k = 1;
% beta_k prior (gamma, rate version)
q_k = 2;
r_k = 0.01;

% Right branch
% mu_dk prior (normal)
f_dk = 0;
g2_dk = 1000;
% tau2_dk prior (inverse gamma)
l_dk = 1;
o_dk = 0.01;

%%% Parameters, will be sampled

% Left branch
xi_ks = n_ks;
upsilon2_ks = zeros(n, 1) + 1000;
alpha_ks = ones(n, 1);
beta_ks = zeros(n, 1) + 0.01;

nu_iks = repmat(xi_ks', N, 1);
sigma2_iks = repmat(upsilon2_ks', N, 1);

% Be careful with p_itks: Subjects have different numbers of trials
trial_nums = histcounts(categorical(subjects(conds == cond))); % This might not be ordered according to subjects
p_itks = zeros(N, max(trial_nums), n);

for i = 1:N
    subject_lens = trial_lens(subjects == i & conds == cond);
    trial_nums(i) = length(subject_lens); % Ensure ordered trial_nums
    
    for t = 1:trial_nums(i)
        p_itks(i, t, :) = initializeBumpIndices(n, subject_lens(t));
    end
end

p_itks = round(p_itks);

% Right branch
mu_dks = zeros(D, n); % Principal component box
tau2_dks = ones(D, n); % Principal component box
m_idks = zeros(N, D, n);

%% Setting up parameter histogram variables

% This block only runs if we want to continue with a previous run. Loads
% the most recent parameter values that were estimated in that run.
if exist('continue_path', 'var') == 1
    load(fullfile(continue_path, 'histograms'), 'hists');
    
    % Which bumps to load
    bump_selection = [1 2 3 6 8];
    
    % Temporary solution for setting start_iter
    start_iter = length(nonzeros(hists.beta_k(:, 1)));
    
    xi_ks = hists.xi_k(start_iter, bump_selection)';
    upsilon2_ks = hists.upsilon2_k(start_iter, bump_selection)';
    alpha_ks = hists.alpha_k(start_iter, bump_selection)';
    beta_ks = hists.beta_k(start_iter, bump_selection)';
    
    nu_iks = squeeze(hists.nu_ik(start_iter, :, bump_selection));
    sigma2_iks = squeeze(hists.sigma2_ik(start_iter, :, bump_selection));
    
    m_idks = squeeze(hists.m_idk(start_iter, :, :, bump_selection));
    mu_dks = squeeze(hists.mu_dk(start_iter, :, bump_selection));
    tau2_dks = squeeze(hists.tau2_dk(start_iter, :, bump_selection));
    
    p_itks = squeeze(hists.p_itk(start_iter, :, :, bump_selection));
    
    % Forward-backward shift to prevent any bump overlap
    p_itks = forwardBackwardShift(p_itks, subjects, conds, trial_lens, ...
        trial_nums, N, n, cond);
end

hists = struct;

hists.xi_k = zeros(max_iter, n);
hists.upsilon2_k = zeros(max_iter, n);
hists.alpha_k = zeros(max_iter, n);
hists.beta_k = zeros(max_iter, n);

hists.nu_ik = zeros(max_iter, N, n);
hists.sigma2_ik = zeros(max_iter, N, n);

hists.m_idk = zeros(max_iter, N, D, n);
hists.mu_dk = zeros(max_iter, D, n);
hists.tau2_dk = zeros(max_iter, D, n);

hists.p_itk = zeros(max_iter, N, max(trial_nums), n);

start_iter = 1;

%% Sampling

tic

for iter = start_iter:max_iter
    disp(iter)
    
    % Record parameter values into histograms variables
    hists.xi_k(iter, :) = xi_ks;
    hists.upsilon2_k(iter, :) = upsilon2_ks;
    hists.alpha_k(iter, :) = alpha_ks;
    hists.beta_k(iter, :) = beta_ks;

    hists.nu_ik(iter, :, :) = nu_iks;
    hists.sigma2_ik(iter, :, :) = sigma2_iks;

    hists.m_idk(iter, :, :, :) = m_idks;
    hists.mu_dk(iter, :, :) = mu_dks;
    hists.tau2_dk(iter, :, :) = tau2_dks;
    
    hists.p_itk(iter, :, :, :) = p_itks;
    
    % Visualize results
    if visualize_flag
        params = struct;
        params.max_iter = max_iter;
        params.iter = iter;
        params.start_iter = start_iter;
        params.normedscore10 = normedscore10;
        params.x = x;
        params.y = y;
        params.subjects = subjects;
        params.coeff10 = coeff10;
        params.latent10 = latent10;
        params.data = data;
        params.chanlocs = chanlocs;
        params.conds = conds;
        params.cond = cond;
        params.relative = false;
        params.firstVis = firstVis;
        
        params.burnIn = 10;
        params.renderLag = 500;
        params.thinLag = 1;

        visualize(params, hists);

        firstVis = false;
    end
    
    disp('Sampling left branch...')
    for k = 1:n
        % For every bump, we sample xi_k, upsilon2_k, alpha_k, and beta_k.
        
        % xi_k
        meanNumerator = s2_k * sum(nu_iks(:, k)) + upsilon2_ks(k) * n_ks(k);
        denominator = N * s2_k + upsilon2_ks(k);
        xi_ks(k) = normrnd(meanNumerator / denominator, ...
            (upsilon2_ks(k) * s2_k) / denominator);
        
        % upsilon2_k
        upsilon2_ks(k) = invgamrnd(N / 2 + a_k, ...
            (sum((nu_iks(:, k) - xi_ks(k)).^2) + 2 * b_k) / 2);
        
        % alpha_k - Pass whole sigma2_iks or just a subset?
        alpha_ks(k) = alpha_kMH(alpha_ks(k), c_k, d_k, e_k, sigma2_iks(:, k), beta_ks(k), 1);
        
        % beta_k
        beta_ks(k) = gamrnd((N + e_k) * a_k + q_k, ... 
            1 / (sum(1 ./ sigma2_iks(:, k)) + r_k));
        
        for i = 1:N
            T = trial_nums(i);
            
            % Then, for every subject, sample nu_ik and sigma2_ik.

            % nu_ik
            meanNumerator = upsilon2_ks(k) * sum(p_itks(i, 1:T, k)) ... 
                + sigma2_iks(i, k) * xi_ks(k);
            denominator = T * upsilon2_ks(k) + sigma2_iks(i, k);
            nu_iks(i, k) = normrnd(meanNumerator / denominator, ...
                (sigma2_iks(i, k) * upsilon2_ks(k)) / denominator);

            % sigma2_ik
            sigma2_iks(i, k) = invgamrnd(T / 2 + alpha_ks(k), ...
                (sum((p_itks(i, 1:T, k) - nu_iks(i, k)).^2) + 2 * beta_ks(k)) / 2);
            
            % Limit sigma2_ik to a minimum - Otherwise this causes crashes
            % in the p_itk sampler.
            sigma2_iks(i, k) = max(sigma2_iks(i, k), 0.1);
            
            % Finally, for every trial, sample p_itk.
            subject_x = x(subjects == i & conds == cond); % Starts of trials
            subject_y = y(subjects == i & conds == cond); % Ends of trials
            
            for t = 1:T
                params = struct;
                params.start = p_itks(i, t, k);
                params.i = i;
                params.t = t;
                params.k = k;
                params.data = normedscore10(subject_x(t):subject_y(t), :);
                params.m_idks = m_idks;
                params.p_itks = p_itks;
                params.nu_iks = nu_iks;
                params.sigma2_iks = sigma2_iks;
                params.P = P;
                params.V = V;
                params.num = 1;
                
                params.mu_dks = mu_dks;
                
                p_itks(i, t, k) = p_itkCalculate(params);
            end
        end
    end
    
    disp('Sampling right branch...')
    for d = 1:D
        for k = 1:n
            % As for the right branch, sample mu_dk and tau^2_dk for every
            % bump and every principal component.
            mu_dks(d, k) = normrnd( ... 
                (g2_dk * sum(m_idks(:, d, k)) + tau2_dks(d, k) * f_dk) / (N * g2_dk + tau2_dks(d, k)), ...
                (tau2_dks(d, k) * g2_dk) / (N * g2_dk + tau2_dks(d, k)) );
            tau2_dks(d, k) = invgamrnd(N / 2 + l_dk, ...
                (sum((m_idks(:, d, k) - mu_dks(d, k)).^2) + 2 * o_dk) / 2);
            
            % Additionally, sample m_idk for every subject.
            for i = 1:N
                % First, acquire all data for this subject.
                subject_x = x(subjects == i & conds == cond);
                
                % T: Number of trials
                T = trial_nums(i);
                
                bump_indices = getBumpIndices(subject_x, p_itks(i, 1:T, k));

                % Get the data
                bumps = normedscore10(bump_indices, d);
                
                bump_corr = repmat(P, 1, T) * bumps;
                
                % And finally draw the sample.
                denominator = 2 * tau2_dks(d, k) * T * sum(P .^ 2) + V;
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
prefix = '';
if dirs{1} == ""
    prefix = '/'; % savePath is absolute
end
parentDir = [prefix fullfile(dirs{1:end-1})];
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

save(fullfile(savePath, 'histograms'), 'hists');

end
