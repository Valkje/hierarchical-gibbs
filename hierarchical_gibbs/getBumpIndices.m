function bump_indices = getBumpIndices(subject_x, p)
%GETBUMPINDICES Calculates all bump indices using p and returns it as a 1D
%array.
%   p contains the middle of all requested bumps, independent of any trial
%   offsets. This function adds those offsets to p and also calculates the
%   bump indices that surround the middle index.
%   - p: Middle index of the bump for every trial, excluding offsets.
%   - subject_x: Offsets for every trial.

bump_indices = p' + subject_x - 1;
bump_indices = repmat(bump_indices, 1, 5) + (-2:2);
bump_indices = reshape(bump_indices', [], 1);

end

