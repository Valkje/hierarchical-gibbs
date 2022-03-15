function bump_indices = getBumpIndices_relative(subject_x, p_itks, i, k)
%GETBUMPINDICES Calculates all bump indices using p and the previous bumps;
%returns it as a 1D array.
%   p_itks contains the middle of all bumps, independent of any trial
%   offsets. To get absolute indices (without offsets), we need to sum the
%   relative values for bumps 1 to k. Then we can add the result to the
%   trial offsets and calculate the bump indices that surround the middle
%   index. 
%   - subject_x: Offsets for every trial.
%   - p_itks: Middle indices of all relative bumps for every trial, 
%     excluding trial offsets. 
%   - i: Subject number.
%   - k: Bump number.

bump_indices = sum(p_itks(i, :, 1:k), 3)' + subject_x - 1;
bump_indices = repmat(bump_indices, 1, 5) + (-2:2);
bump_indices = reshape(bump_indices', [], 1);

end

