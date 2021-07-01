function samples = alpha_kMH(start, a_k, b_k, q_k, alpha_iks, beta_k, n)
%ALPHA_KMH Summary of this function goes here
%   Detailed explanation goes here
    % 10 samples burn-in, 1 sample autocorrelation lag
    burnIn = 10;
    lag = 10;
    skipCounter = burnIn + lag;
    
    % Starting value of 2
    xPrevious = start;
    samples = zeros(n, 1);
    acceptedNum = 0;
    
    % Variance of proposal distribution
    proposal_variance = 2;
    
    product = prod(alpha_iks);

    while true
        xProposed = normrnd(xPrevious, 2, 1);
        
        % Constraints on xProposed
        if xProposed <= 0
            continue % Try again
        end
        
        acceptance = alpha_kPdf(xProposed, a_k, b_k, q_k, alpha_iks, beta_k, product) ... 
            / alpha_kPdf(xPrevious, a_k, b_k, q_k, alpha_iks, beta_k, product);
        % Correcting for the truncated proposal distribution
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

