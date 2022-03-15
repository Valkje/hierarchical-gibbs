function samples = invgamrnd(alpha, beta, n)
%INVGAMRND Generates values from an inverse gamma distribution.
%   Uses Matlab's gamma distribution sampler to calculate a sample for the
%   inverse gamma distribution. Note that beta corresponds to the rate
%   parameter of the gamma distribution, while the Matlab function requires
%   a scale parameter. We fix this discrepancy by simply calculating 1 /
%   beta.
arguments
   alpha double
   beta double
   n double = 1
end

    samples = 1 ./ gamrnd(alpha, 1 / beta, n, 1);
end

