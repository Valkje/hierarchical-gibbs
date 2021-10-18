function x = dirrnd(alpha, n)
%DIRRND Samples from a Dirichlet distribution
%   Uses gamma samplers to sample from a Dirichlet distribution.
%   - alpha: Parameter vector for the Dirichlet distribution.
%   - n: Number of samples to draw.
    k = length(alpha);
    gammas = gamrnd(repmat(alpha, n, 1), 1, n, k);
    x = gammas ./ repmat(sum(gammas, 2), 1, k);
end

