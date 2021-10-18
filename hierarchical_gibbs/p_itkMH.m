function samples = p_itkMH(params)
%P_ITKMH Creates num samples for p_itk.
%   A Metropolis sampler for p_itk.

    start = params.start;
    i = params.i;
    t = params.t;
    k = params.k;
    data = params.data';
    m_idks = params.m_idks;
    p_itks = params.p_itks;
    nu_iks = params.nu_iks;
    sigma2_iks = params.sigma2_iks;
    P = params.P;
    V = params.V;
    num = params.num;

    burnIn = 10;
    lag = 0;
    skipCounter = burnIn + lag;
    
    % Set starting value
    xPrevious = start;
    samples = zeros(num, 1);
    numAccepted = 0;
    
    % Variance of proposal distribution
    proposal_variance = 10;

    while true
        xProposed = round(normrnd(xPrevious, proposal_variance, 1));
        
        % Constraints on xProposed - truncated proposal distribution
        if xProposed <= 0
            continue % Try again
        end
        
        acceptance = p_itkPdf2(xProposed,  i, t, k, data, m_idks, p_itks, nu_iks, sigma2_iks, P, V) ... 
            / p_itkPdf2(xPrevious,  i, t, k, data, m_idks, p_itks, nu_iks, sigma2_iks, P, V);
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
            numAccepted = numAccepted + 1;
            samples(numAccepted) = xPrevious;
        end

        if numAccepted == num
            break
        end
    end

end

