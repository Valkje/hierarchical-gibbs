function indices = initializeBumpIndices(n, trial_len)
%INITIALIZEBUMPINDICES Summary of this function goes here
%   Detailed explanation goes here

flat_proportions = dirrnd(2, n+1);

indices = round(trial_len * cumsum(flat_proportions(1:end-1)));

% Forward shift
if indices(1) < 3
    indices(1) = 3;
end

for i = 2:n
    if indices(i) < indices(i-1) + 5
        indices(i) = indices(i-1) + 5;
    end
end

% Backward shift
if indices(end) > trial_len - 3
    indices(end) = trial_len - 3;
end

for i = n-1:-1:1
    if indices(i) > indices(i+1) - 5
        indices(i) = indices(i+1) - 5;
    end
end

end

