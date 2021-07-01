function samples = tiMH(data, t_is, alphas, V, x, y, numBumps, n)
%TIMH Summary of this function goes here
%   Detailed explanation goes here
    
    % Sum all the PCs together
    data = sum(data, 2);
    lens = y - x + 1;
    
    % Calculate alpha_tildes
    alpha_tildes = zeros(length(alphas), 1);
    for i = 1:length(alphas)
        alpha_tildes(i) = sum(alphas(i+1:end));
    end
    
    % Convert t_is to z_is
    ys = sqrt(t_is);

    angles = zeros(1, length(ys) - 1);
    for i = 1:length(angles)
        angles(i) = acos(ys(i) / sqrt(sum( ys(i:end).^2 )));
    end

    z_is = sin(angles).^2;
    
    % Minimal duration is equivalent to 1 sample
    min_t = 1 / min(lens);
    
    % Chain setup
    samples = zeros(n, length(t_is));

    burnin = 50;
    lag = 10;
    skipCounter = burnin + lag;

    i = 0;
    k = 1;

    % Note: We are not sampling t_is(end), as it is fully determined by the
    % other t_i's

    while true
        % Calculate minimal z_k
%         t_is_temp = t_is;
%         t_is_temp(k) = min_t;
%         
%         ys = sqrt(t_is_temp);
%         angle = acos(ys(k) / sqrt(sum( ys(k:end).^2 )));
%         max_z = sin(angle).^2;
        max_z = 1;
        
        % Propose new z
        zProposed = unifrnd(0, max_z, 1);

        % Calculate proposed t_is
        z_is_prop = z_is;
        z_is_prop(k) = zProposed;
        t_is_prop = zeros(1, numBumps+1);
        for j = 1:length(alphas)-1
            t_is_prop(j) = prod(z_is_prop(1:j-1)) * (1 - z_is_prop(j));
        end
        t_is_prop(end) = prod(z_is_prop);
        
        
        [starts, ends] = getFlatIndices(x, y, numBumps, t_is);
        [starts_prop, ends_prop] = ...
            getFlatIndices(x, y, numBumps, t_is_prop);

        S_sqr_sum = 0;
        S_sqr_sum_prop = 0;
        for h = 1:length(x)
            S_sqr_sum = S_sqr_sum + sum( data(starts(h):ends(h)).^2 );
            S_sqr_sum_prop = S_sqr_sum_prop + sum( data(starts_prop(h):ends_prop(h)).^2 );
        end

        % Avoid NaNs
        log_numerator = log(zProposed ^ (alpha_tildes(k) - 1) * (1 - zProposed) ^ (alphas(k) - 1)) + -S_sqr_sum_prop / V;
        log_denominator = log(z_is(k) ^ (alpha_tildes(k) - 1) * (1 - z_is(k)) ^ (alphas(k) - 1)) + -S_sqr_sum / V;
        acceptance = exp(log_numerator - log_denominator);
        
%         acceptance = zProposed ^ (alpha_tildes(k) - 1) * (1 - zProposed) ^ (alphas(k) - 1) * exp(-S_sqr_sum_prop / V) ... 
%             / (z_is(k) ^ (alpha_tildes(k) - 1) * (1 - z_is(k)) ^ (alphas(k) - 1) * exp(-S_sqr_sum / V));

        if unifrnd(0, 1, 1) <= acceptance
            z_is(k) = zProposed;
            t_is = t_is_prop;
        end

        k = k + 1;
        if k == length(t_is)
            t_is(k) = 1 - sum(t_is(1:end-1));

            k = 1;

            if skipCounter > 0
                skipCounter = skipCounter - 1;
            else
                skipCounter = lag;
                i = i + 1;
                samples(i, :) = t_is;
            end
        end

        if i == n
            break
        end
    end
end

