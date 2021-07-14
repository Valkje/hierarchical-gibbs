function samples = alpha_kMH(start, a_k, b_k, q_k, alpha_iks, beta_k, n)
%ALPHA_KMH Creates n samples for alpha_k.
%   A Metropolis sampler for alpha_k. Has a burn-in period of 10 samples
%   and and additional lag of 10 samples to thin the chain. Uses a normal
%   proposal distribution with variance 2, centered on the previous value
%   in the chain.


    % 10 samples burn-in, 10 samples autocorrelation lag
    burnIn = 10;
    lag = 10;
    skipCounter = burnIn + lag;
    
    xPrevious = start;
    samples = zeros(n, 1);
    acceptedNum = 0;
    
    % Variance of proposal distribution
    proposal_variance = 2;
    
    product = prod(alpha_iks);

    while true
        xProposed = normrnd(xPrevious, proposal_variance, 1);
        
        % Constraints on xProposed
        if xProposed <= 0
            continue % Try again
        end
        
        acceptance = alpha_kPdf(xProposed, a_k, b_k, q_k, alpha_iks, beta_k, product) ... 
            / alpha_kPdf(xPrevious, a_k, b_k, q_k, alpha_iks, beta_k, product);
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
            acceptedNum = acceptedNum + 1;
            samples(acceptedNum) = xPrevious;
        end

        if acceptedNum == n
            break
        end
    end
end

