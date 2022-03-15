function indices = initializeBumpIndices(n, trial_len)
%INITIALIZEBUMPINDICES Initializes bump indices using a Dirichlet
%distribution.
%   Samples a vector of size n whose elements sum to 1, then scales every
%   element by the length of the trial. For more explanation of the forward
%   and backward shifting, please see FORWARDBACKWARDSHIFT.

flat_proportions = dirrnd(ones(1, n+1), 1);

indices = round(trial_len * cumsum(flat_proportions(1:end-1)));

%% This section should become part of FORWARDBACKWARDSHIFT
% Forward shift
if indices(1) < 3
    indices(1) = 3;
end

for i = 2:n
    if indices(i) < indices(i-1) + 1
        indices(i) = indices(i-1) + 1;
    end
end

% Backward shift
if indices(end) > trial_len - 3
    indices(end) = trial_len - 3;
end

for i = n-1:-1:1
    if indices(i) > indices(i+1) - 1
        indices(i) = indices(i+1) - 1;
    end
end

end

