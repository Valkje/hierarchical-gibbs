function density = p_itkPdf2(p_itk, i, t, k, data, m_idks, p_itks, nu_iks, sigma2_iks, P, V)
%p_itkPdf Provides an alternative pdf of $p_{itk}$.
%   Provides the probability density function of $p_{itk}$. Requires DATA
%   to consist of only the D-dimensional data for participant i and trial
%   t.

n = size(p_itks, 3);
trial_length = size(data, 1);

% Account for possible illegal values
if p_itk < 3 || p_itk > trial_length - 2
    density = 0;
    return
end

flat_indices = true(trial_length, 1);

% Mark all bump indices as non-flat
for p = 1:n
    if p == k
        flat_indices(p_itk-2 : p_itk+2) = false;
    % Check that the bumps do not overlap
    elseif p_itk > p_itks(i, t, p) - 5 && p_itk < p_itks(i, t, p) + 5
        density = 0;
        return
    end
    
    flat_indices(p_itks(i, t, p)-2 : p_itks(i, t, p)+2) = false;
end

flat_sum = sum( data(flat_indices, :).^2, 'all' );

% Bump itself
% bump_sum = sum( data(p_itk-2:p_itk+2, :).^2 - ... 
%     2 * data(p_itk-2:p_itk+2, :) .* m_idks(i, :, k) .* P', 'all' );

% Experimentation with simple correlation
bump_sum = -sum( data(p_itk-2:p_itk+2, :) .* m_idks(i, :, k) .* P', 'all' );

density = exp(- (flat_sum + bump_sum) / V) * exp(-(p_itk^2 - 2*p_itk*nu_iks(i, k) + nu_iks(i, k)^2) / (2 * sigma2_iks(i, k)));
end

