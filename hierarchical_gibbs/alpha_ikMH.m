function samples = alpha_ikMH(alpha_iks, start, t_ik, alpha_k, beta_k, n)
%ALPHA_IKMH Creates n samples for alpha_ik.
%   A Metropolis sampler for alpha_ik. Has a burn-in period of 20 samples
%   and and additional lag of 10 samples to thin the chain. Uses a normal
%   proposal distribution with variance 2, centered on the previous value
%   in the chain.
    % 20 samples burn-in, 10 samples autocorrelation lag
    burnIn = 20;
    lag = 10;
    skipCounter = burnIn + lag;
    
    % Starting value of 2
    xPrevious = start;
    samples = zeros(n, 1);
    i = 0;
    
    % Variance of proposal distribution
    proposal_variance = 2;

    while true
        xProposed = normrnd(xPrevious, proposal_variance, 1);
        
        % Constraints on xProposed - truncated proposal distribution
        if xProposed <= 0
            continue % Try again
        end
        
        acceptance = alpha_ikPdf(xProposed, alpha_iks, t_ik, alpha_k, beta_k) ... 
            / alpha_ikPdf(xPrevious, alpha_iks, t_ik, alpha_k, beta_k);
        % Correcting for the truncated proposal distribution
        % See https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/
        acceptance = acceptance * normcdf(xPrevious/sqrt(proposal_variance)) ...
            / normcdf(xProposed/sqrt(proposal_variance));
        
        if unifrnd(0, 1, 1) <= acceptance
            xPrevious = xProposed;
        end
        
        if skipCounter > 0
            skipCounter = skipCounter - 1;
        else
            skipCounter = lag;
            i = i + 1;
            samples(i) = xPrevious;
        end

        if i == n
            break
        end
    end
end

