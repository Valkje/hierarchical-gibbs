function [] = main_bumps(path, save_dir, dataset, visualize_flag, parameters)
% Gibbs sampler for discovering cognitive stages in EEG data

% There are many variables that need to be sampled in turn:
% - We start the entire procedure by fixing every variable to some non-
% specific value. (Starting values should, in principle, not matter.)
% - Then, repeat until convergence:
%   - Sample every variable within the left branch. For every bump, we sample
%   xi_k, upsilon2_k, alpha_k, and beta_k. Then, for every subject, sample
%   nu_ik and sigma2_ik. Finally, for every trial, sample p_itk. 
%   - As for the right branch, sample mu_dk and tau^2_dk for every
%   bump and every principal component. Additionally, sample m_idk for
%   every subject.

arguments
    path char = '~/Documents/Uni/FYRP/'
    save_dir char = 'data/run_debug'
    dataset string = "assoc_recog_overall_SNR5"
    visualize_flag logical = false
    parameters struct = struct
end

%% Checking validity of and setting defaults for PARAMETERS

if ~isfield(parameters, 'max_bumps')
    parameters.max_bumps = 10; % Number of bumps
end

if ~isfield(parameters, 'N')
    parameters.N = 2; % Number of subjects
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

% Whether to continue with a previous run
% continue_path = [path 'data/run2_debug_4'];

%% Check what happens when we delete some bumps

% Delete bump 4
% delete_indices = all_pos(:, 4) + x - 1;
% delete_indices = reshape(delete_indices, [], 1);
% delete_indices = repmat(delete_indices, 1, 5) + (-2:2);
% delete_indices = reshape(delete_indices, [], 1);
% normedscore10(delete_indices, :) = normrnd(0, 1, length(delete_indices), 10);

%% Setting up EEGlab

if visualize_flag
    eeglab;
    close % the EEGlab figure
end

%% Inverting every trial

% for trl = 1:length(x)
%     normedscore10(x(trl):y(trl), :) = normedscore10(y(trl):-1:x(trl), :);
% end

%% Cutting off the first two bumps - requires adjusting x, y, and trial_lengths

% newX = zeros(length(x), 1);
% newY = zeros(length(y), 1);
% for trl = 1:length(x)
%     trimmed = normedscore10(x(trl) + all_pos(trl, 2) + 2:y(trl), :);
%     
%     if trl == 1
%         newX(1) = 1;
%         newY(1) = length(trimmed);
%     else
%         newX(trl) = newY(trl-1) + 1;
%         newY(trl) = newX(trl) + length(trimmed) - 1;
%     end
%     
%     newNormedscore10(newX(trl):newY(trl), :) = trimmed;
% end
% 
% x = newX;
% y = newY;
% normedscore10 = newNormedscore10;
% 
% trial_lens = y - x + 1;

%% What happens if we make bump 3 negative?

% Bump 3 positions
% bump3Pos = x + all_pos(:, 3) - 1;
% bump3Pos = bump3Pos + repmat(-2:2, length(bump3Pos), 1);
% bump3Pos = reshape(bump3Pos', [], 1);
% 
% normedscore10(bump3Pos, :) = -normedscore10(bump3Pos, :);

%% Setting up variables

max_iter = 10;

% Used to initialize figures
firstVis = true;

% Which condition to sample for.
cond = 1;

% A normalising variable from Anderson et al. (2016)
V = 5;

% N = length(unique(subjects));
N = parameters.N;
D = size(normedscore10, 2);

% Bump topology - 1 by 5 matrix
P = sin(pi/10:pi/5:pi-pi/10);

max_bumps = parameters.max_bumps;

% xi_k prior (normal)
s2_k = 100;
% upsilon2_k prior (inverse gamma)
a_k = 1;
b_k = 0.01;
% alpha_k prior (non-closed form)
c_k = 50; % was 50, 2021-11-01
d_k = 0.5;
e_k = 1;
% beta_k prior (gamma, rate version)
q_k = 2;
r_k = 0.01; % was 0.001, 2021-11-01

% Right branch
% mu_dk prior (normal)
f_dk = 0;
g2_dk = 1000;
% tau2_dk prior (inverse gamma)
l_dk = 1;
o_dk = 0.01;


% Matrices for transitioning between different values of n
n_ks_max = ((1:max_bumps) * (mean(trial_lens) / (max_bumps + 1)))';
xi_ks_max = n_ks_max;
upsilon2_ks_max = zeros(max_bumps, 1) + 1000;
alpha_ks_max = ones(max_bumps, 1);
beta_ks_max = zeros(max_bumps, 1) + 0.01;

nu_iks_max = repmat(xi_ks_max', N, 1);
sigma2_iks_max = repmat(upsilon2_ks_max', N, 1);

% Be careful with p_itks: Subjects have different numbers of trials
trial_nums = histcounts(categorical(subjects(conds == cond))); 
p_itks_max = zeros(N, max(trial_nums), max_bumps);

% Right branch
mu_dks_max = zeros(D, max_bumps); % Principal component box
tau2_dks_max = ones(D, max_bumps); % Principal component box
m_idks_max = zeros(N, D, max_bumps);

% Histograms
hists = struct;

hists.xi_k = zeros(max_iter, max_bumps);
hists.upsilon2_k = zeros(max_iter, max_bumps);
hists.alpha_k = zeros(max_iter, max_bumps);
hists.beta_k = zeros(max_iter, max_bumps);

hists.nu_ik = zeros(max_iter, N, max_bumps);
hists.sigma2_ik = zeros(max_iter, N, max_bumps);

hists.m_idk = zeros(max_iter, N, D, max_bumps);
hists.mu_dk = zeros(max_iter, D, max_bumps);
hists.tau2_dk = zeros(max_iter, D, max_bumps);

hists.p_itk = zeros(max_iter, N, max(trial_nums), max_bumps);

for n = 1:max_bumps
    
    % Left branch
    % xi_k prior (normal)
    n_ks = ((1:n) * (mean(trial_lens) / (n + 1)))'; % Even spread of bumps as prior

    %%% Parameters, will be sampled

    % Left branch
    xi_ks = n_ks;
    upsilon2_ks = zeros(n, 1) + 1000;
    alpha_ks = ones(n, 1);
    beta_ks = zeros(n, 1) + 0.01;

    nu_iks = repmat(xi_ks', N, 1);
    sigma2_iks = repmat(upsilon2_ks', N, 1);

    % Be careful with p_itks: Subjects have different numbers of trials
    trial_nums = histcounts(categorical(subjects(conds == cond))); 
    p_itks = zeros(N, max(trial_nums), n);

    for i = 1:N
        subject_lens = trial_lens(subjects == i & conds == cond);
        trial_nums(i) = length(subject_lens); % Ensure ordered trial_nums

        for t = 1:trial_nums(i)
            % All p_itks for a given k are relatively close right now; perhaps
            % not a fair initialisation
            p_itks(i, t, :) = (1:n) * (subject_lens(t) / (n + 1));
        end
    end

    p_itks = round(p_itks);

    % Right branch
    mu_dks = zeros(D, n); % Principal component box
    tau2_dks = ones(D, n); % Principal component box
    m_idks = zeros(N, D, n);
    
    
    new_bump_k = 1;
    
    if n ~= 1
        % Determine most optimal spot for next bump:
        % - loop over all flat proportions
        % - select flat based on largest average proportion
        flat_props = zeros(N, n);
        for k = 1:n % k stands for flat
            for i = 1:N
                % Note: subject_lens is a column vector
                subject_lens = trial_lens(subjects == i & conds == cond);
                
                % How we calculate the proportions of a flat depends on
                % whether the flat is the first one (k == 1), the last one
                % (k == n), or any flat in-between.
                if k == 1
                    % The location of the first bump defines the length of
                    % the first flat.
                    subject_ps = ...
                        squeeze(p_itks_max(i, 1:trial_nums(i), k))';
                    flat_props(i, k) = mean(subject_ps ./ subject_lens);
                elseif k == n
                    % The length of the last flat is defined by the length
                    % of the trial minus the position of the last bump.
                    prev_subject_ps = ...
                        squeeze(p_itks_max(i, 1:trial_nums(i), k-1))';
                    flat_props(i, k) = ...
                        mean((subject_lens - prev_subject_ps) ./ subject_lens);
                else
                    % For all the other flats, the length is defined by the
                    % location of the current bump minus the location of
                    % the previous bump.
                    subject_ps = ...
                        squeeze(p_itks_max(i, 1:trial_nums(i), k))';
                    prev_subject_ps = ...
                        squeeze(p_itks_max(i, 1:trial_nums(i), k-1))';
                    flat_props(i, k) = ...
                        mean((subject_ps - prev_subject_ps) ./ subject_lens);
                end
            end
        end
        
        % Find the flat that is, on average, the longest. The number of
        % that flat equals the number of the bump we will insert.
        [~, new_bump_k] = max(mean(flat_props, 1));
        
        % Although we have determined the flat that has the most space on
        % average, in an individual trial that flat might be too short to
        % harbor any new bumps. Therefore, we first need to calculate the
        % available room for the new bump for every individual trial. Next,
        % if it turns out there is too little room available, we have to
        % shift the already existing bumps to make that room.
        
        new_locations = zeros(N, max(trial_nums));
        
        for i = 1:N
            subject_lens = trial_lens(subjects == i & conds == cond);
            
            for t = 1:trial_nums(i)
                indices = p_itks_max(i, t, :); % Just a shorter alias
                trial_len = subject_lens(t);
                
                if new_bump_k == 1
                    % Difference between following bump and start
                    room_for_bump = indices(new_bump_k) - 3;
                    % While we're at it, also calculate the offset for new
                    % bump placement
                    offset = 2;
                elseif new_bump_k == n
                    % Difference between preceding bump and end
                    room_for_bump = trial_len - indices(new_bump_k-1) - 2;
                    offset = indices(new_bump_k - 1);
                else
                    % Difference between two bumps
                    room_for_bump = indices(new_bump_k) - ...
                        indices(new_bump_k-1) - 1;
                    offset = indices(new_bump_k - 1);
                end
                
                required_room = 1;
                if room_for_bump < required_room
                    % Make space for the new bump
                    if new_bump_k == n
                        % Shift left bump before new bump
                        p_itks_max(i, t, new_bump_k-1) = ...
                            indices(new_bump_k-1) - (required_room - room_for_bump);
                    else
                        % Shift right bump after new bump
                        p_itks_max(i, t, new_bump_k) = ...
                            indices(new_bump_k) + (required_room - room_for_bump);
                    end
                    
                    % Any forward-backward shifting due to bumps that might
                    % now overlap should occur after the new bump has
                    % actually been placed
                    
                    % If we're here, we have made some room
                    room_for_bump = required_room;
                end
                
                % We will now randomly place the new bump within the space
                % that is available for it, given the calculated offset and
                % room.
                new_locations(i, t) = offset + randi(room_for_bump);
            end
        end
        
        % Load all the 'old' bumps from the transitioning matrices
        load_sel = true(n, 1);
        load_sel(new_bump_k) = false;
        
        xi_ks(load_sel) = xi_ks_max(1:n-1);
        upsilon2_ks(load_sel) = upsilon2_ks_max(1:n-1);
        alpha_ks(load_sel) = alpha_ks_max(1:n-1);
        beta_ks(load_sel) = beta_ks_max(1:n-1);
        
        nu_iks(:, load_sel) = nu_iks_max(:, 1:n-1);
        sigma2_iks(:, load_sel) = sigma2_iks_max(:, 1:n-1);
        
        p_itks(:, :, load_sel) = p_itks_max(:, :, 1:n-1);
        
        mu_dks(:, load_sel) = mu_dks_max(:, 1:n-1);
        tau2_dks(:, load_sel) = tau2_dks_max(:, 1:n-1);
        m_idks(:, :, load_sel) = m_idks_max(:, :, 1:n-1);
        
        % Make space within the hists field for the new bump values
        hists.xi_k(:, load_sel) = hists.xi_k(:, 1:n-1);
        hists.upsilon2_k(:, load_sel) = hists.upsilon2_k(:, 1:n-1);
        hists.alpha_k(:, load_sel) = hists.alpha_k(:, 1:n-1);
        hists.beta_k(:, load_sel) = hists.beta_k(:, 1:n-1);

        hists.nu_ik(:, :, load_sel) = hists.nu_ik(:, :, 1:n-1);
        hists.sigma2_ik(:, :, load_sel) = hists.sigma2_ik(:, :, 1:n-1);

        hists.m_idk(:, :, :, load_sel) = hists.m_idk(:, :, :, 1:n-1);
        hists.mu_dk(:, :, load_sel) = hists.mu_dk(:, :, 1:n-1);
        hists.tau2_dk(:, :, load_sel) = hists.tau2_dk(:, :, 1:n-1);
        
        hists.p_itk(:, :, :, load_sel) = hists.p_itk(:, :, :, 1:n-1);
        
        % Actually place the new bumps
        for i = 1:N
            p_itks(i, 1:trial_nums(i), new_bump_k) = ...
                new_locations(i, 1:trial_nums(i));
        end
        
        % We might have created some overlapping bumps by now, so let's
        % make sure they do not longer overlap by forward-backward shifting
        p_itks = forwardBackwardShift(p_itks, subjects, conds, ...
            trial_lens, trial_nums, N, n, cond);
    end
    
    % Sampling
    start_iter = 1;
    for iter = start_iter:max_iter
        disp(iter)

        % Record parameter values into histograms variables
        hists.xi_k(iter, new_bump_k) = xi_ks(new_bump_k);
        hists.upsilon2_k(iter, new_bump_k) = upsilon2_ks(new_bump_k);
        hists.alpha_k(iter, new_bump_k) = alpha_ks(new_bump_k);
        hists.beta_k(iter, new_bump_k) = beta_ks(new_bump_k);

        hists.nu_ik(iter, :, new_bump_k) = nu_iks(:, new_bump_k);
        hists.sigma2_ik(iter, :, new_bump_k) = sigma2_iks(:, new_bump_k);

        hists.m_idk(iter, :, :, new_bump_k) = m_idks(:, :, new_bump_k);
        hists.mu_dk(iter, :, new_bump_k) = mu_dks(:, new_bump_k);
        hists.tau2_dk(iter, :, new_bump_k) = tau2_dks(:, new_bump_k);

        hists.p_itk(iter, :, :, new_bump_k) = p_itks(:, :, new_bump_k);

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
            params.renderLag = 100;
            params.thinLag = 1;
            
            params.n = n;
            
            visualize(params, hists);
            
            firstVis = false;
        end

        disp('Sampling left branch...')
        for k = new_bump_k
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
                
                % Limit sigma2_ik to a minimum
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
        
        % Setting bump positions to their mode before sampling the final
        % topologies
        if iter == max_iter
            for i = 1:N
                for t = 1:trial_nums(i)
                    [counts, edges] = histcounts(hists.p_itk(:, i, t, new_bump_k), ...
                        range(hists.p_itk(:, i, t, new_bump_k)) + 1);
                    [~, modeIdx] = max(counts);
                    p_itks(i, t, new_bump_k) = round(...
                        mean([edges(modeIdx) edges(modeIdx+1)]));
                end
            end
        end

        disp('Sampling right branch...')
        for d = 1:D
            for k = new_bump_k
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

                    % Get those data!
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
    
    % Save to transitioning matrices
    xi_ks_max(1:n) = xi_ks;
    upsilon2_ks_max(1:n) = upsilon2_ks;
    alpha_ks_max(1:n) = alpha_ks;
    beta_ks_max(1:n) = beta_ks;
    
    nu_iks_max(:, 1:n) = nu_iks;
    sigma2_iks_max(:, 1:n) = sigma2_iks;
    
    % Experiment with saving the mode instead of the value that was sampled
    % last
    for i = 1:N
        for t = 1:trial_nums(i)
            [counts, edges] = histcounts(hists.p_itk(:, i, t, new_bump_k), ...
                range(hists.p_itk(:, i, t, new_bump_k)) + 1);
            [~, modeIdx] = max(counts);
            p_itks(i, t, new_bump_k) = round(...
                mean([edges(modeIdx) edges(modeIdx+1)]));
        end
    end
    p_itks_max(:, :, 1:n) = p_itks;
    
    mu_dks_max(:, 1:n) = mu_dks;
    tau2_dks_max(:, 1:n) = tau2_dks;
    m_idks_max(:, :, 1:n) = m_idks;
end

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
