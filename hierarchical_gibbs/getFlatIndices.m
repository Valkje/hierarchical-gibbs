function [starts_offset, ends_offset] = getFlatIndices(x, y, n, t_is)
%GETFLATINDICES Summary of this function goes here
%   Detailed explanation goes here

    lens = y - x + 1;

    % C: Number of trials
    C = length(x);

    % Subtract the number of samples reserved for the bumps
    total_flat_lens = lens - n * 5;

%     % Then, for every trial, get every concrete flat duration
%     % by multiplying proportions with the number of samples
%     % reserved for flats for every trial.
%     flat_lens = repmat(total_flat_lens, 1, n + 1) .* ...
%         repmat(t_is, C, 1);
%     flat_lens = round(flat_lens);
%     
%     % We will not allow flats with length 0 - give them length 1, and
%     % subtract the necessary amount from the final flat.
%     zero_indices = flat_lens == 0;
%     flat_lens(zero_indices) = 1;
%     flat_lens(:, end) = flat_lens(:, end) - sum(zero_indices, 2);
% 
%     flat_starts = ones(size(flat_lens));
%     flat_starts(:, 2:end) = flat_lens(:, 1:end-1) + 5;
%     flat_starts = cumsum(flat_starts, 2);
% 
%     starts_offset = flat_starts + x - 1;
% 
%     ends_offset = starts_offset + flat_lens - 1;
% 
%     % Due to rounding errors, the end index of every trial might not
%     % match the index found in y. Therefore, we truncate or extend all
%     % end indices
%     ends_offset(:, end) = y;
    
    %% Temporary - bumps at fixed locations
    
%     min_flat_len = 77 - n * 5;
    min_flat_len = min(total_flat_lens);
    flat_lens = repmat(min_flat_len, 1, n + 1) .* t_is;
    flat_lens = round(flat_lens);
    
    zero_indices = flat_lens == 0;
    flat_lens(zero_indices) = 1;
    flat_lens(:, end) = flat_lens(:, end) - sum(zero_indices, 2);
    
    flat_starts = ones(size(flat_lens));
    flat_starts(:, 2:end) = flat_lens(:, 1:end-1) + 5;
    flat_starts = cumsum(flat_starts, 2);
    flat_starts = repmat(flat_starts, C, 1);
    
    starts_offset = flat_starts + x - 1;

    ends_offset = starts_offset + flat_lens - 1;
    
    ends_offset(:, end) = y;
end

