function bump_indices = getBumpIndices(subject_x, p)
%GETBUMPINDICES Summary of this function goes here
%   p - Middle index of the bump for every trial. Excludes offset.
%   subject_x - Offsets for every trial.

bump_indices = p' + subject_x - 1;
bump_indices = repmat(bump_indices, 1, 5) + (-2:2);
bump_indices = reshape(bump_indices', [], 1);

%     % Subtract the number of samples reserved for the bumps
%     total_flat_lens = subject_lens - n * 5;
% 
%     % Then, for every trial, get every concrete flat duration
%     % by multiplying proportions with the number of samples
%     % reserved for flats for every trial.
%     flat_lens = repmat(total_flat_lens, 1, n + 1) .* ...
%         repmat(ts, length(total_flat_lens), 1);
%     flat_lens = round(flat_lens);
% 
%     % Now calculate the cumulative sum - after re-adding the
%     % number of preceding bump samples - to get the sample that
%     % signifies the end of every flat for every trial.
%     flat_ends = flat_lens;
%     flat_ends(:, 2:n+1) = flat_ends(:, 2:n+1) + 5;
%     flat_ends = cumsum(flat_ends, 2);
% 
%     % Finally, determine the bump locations for all trials for
%     % flat k.
%     bump_starts = flat_ends(:, k) + 1;
%     
%     % To index into the entire data array, we need to offset
%     % the bump locations by every trial's onset.
%     bump_starts = subject_x + bump_starts - 1; % C by 1
% 
%     % Determine bump indices, reshape into column vector.
%     bump_indices = repmat(bump_starts, 1, 5) + (0:4);
%     bump_indices = reshape(bump_indices', [], 1);
    
    %% Fixed bump positions
    
%     C = length(subject_x);
%     
%     min_flat_len = min(total_flat_lens);
%     flat_lens = repmat(min_flat_len, 1, n + 1) .* ts;
%     flat_lens = round(flat_lens);
%     
%     flat_ends = flat_lens;
%     flat_ends(:, 2:n+1) = flat_ends(:, 2:n+1) + 5;
%     flat_ends = cumsum(flat_ends, 2);
%     
%     bump_starts = flat_ends(:, k) + 1;
%     bump_starts = repmat(bump_starts, C, 1);
%     
%     bump_starts = subject_x + bump_starts - 1;
%     
%     bump_indices = repmat(bump_starts, 1, 5) + (0:4);
%     bump_indices = reshape(bump_indices', [], 1);
end

