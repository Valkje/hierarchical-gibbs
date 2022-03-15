function samples = p_itkCalculate_relative(params)
%P_ITKCALCULATE_RELATIVE Calculates a probability per trial index and
%samples indices accordingly.
%   Calculates a density for every index of the trial (using the relative
%   pdf), and then normalizes all those densities to get probabilities.
%   Next, samples are drawn from a uniform distribution on (0, 1) and
%   interp1 is used to assign those numbers to trial indices according to
%   the probabilities that have been calculated.

% T: Trial size
T = size(params.data, 1);

densities = zeros(T, 1);
for idx = 1:T
    densities(idx) = p_itkPdf_relative(idx, params);
end

probabilities = densities ./ sum(densities);

% It can be that nu is too distant from any viable position - then
% probabilities will consist of NaNs. In that case, just return the
% starting value.
if all(isnan(probabilities))
    samples = repmat(params.start, params.num, 1);
    %disp('Ack!')
    return
end

% See https://stackoverflow.com/questions/31665504/generate-random-samples-from-arbitrary-discrete-probability-density-function-in/
cdf = cumsum(probabilities); %remember, first term out of of cumsum is not zero

% Because of the operation we're doing below (interp1 followed by ceil)
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

