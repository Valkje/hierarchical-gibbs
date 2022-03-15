function samples = p_itkCalculate(params)
%P_ITKCALCULATE Calculates a probability per trial index and samples
%indices accordingly.
%   Calculates a density for every index of the trial, and then normalizes
%   all those densities to get probabilities. Next, samples are drawn from
%   a uniform distribution on (0, 1) and interp1 is used to assign those
%   numbers to trial indices according to the probabilities that have been
%   calculated.

% T: Trial size
T = size(params.data, 1);

densities = zeros(T, 1);
for idx = 1:T
    % p_itkPdf: Standard
    % p_itkPdf2: Jumping bumps
    densities(idx) = p_itkPdf(idx, params);
end

probabilities = densities ./ sum(densities);

% See https://stackoverflow.com/questions/31665504/generate-random-samples-from-arbitrary-discrete-probability-density-function-in/
cdf = cumsum(probabilities); %remember, first term out of of cumsum is not zero

% Because of the operation we're doing below (interp1)
% we need the CDF to start at zero
cdf = [0; cdf(:)];

% Generate random values
rand_vals = rand(params.num, 1);  %spans zero to one

% Look into CDF to see which index the rand val corresponds to
% Take into account that some probabilities might be zero
[~, uniqInd] = unique(cdf);
indexArr = 0:1/(length(cdf)-1):1;
out_val = interp1(cdf(uniqInd), indexArr(uniqInd), rand_vals, 'next');

samples = round(out_val * T);

end

