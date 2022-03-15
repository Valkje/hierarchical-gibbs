function density = p_itkPdf(p_itk, params)
%p_itkPdf Provides the pdf of p_itk.
%   Provides the probability density function of p_itk. Requires data
%   to consist of only the D-dimensional data for participant i and trial
%   t.

i = params.i;
t = params.t;
k = params.k;
data = params.data;
m_idks = params.m_idks;
p_itks = params.p_itks;
nu_iks = params.nu_iks;
sigma2_iks = params.sigma2_iks;
P = params.P;
V = params.V;
mu_dks = params.mu_dks;

% Hackedyhackhack: Invert the sign of m_idk if it differs from mu_dk
m_idks(i, :, k) = sign(m_idks(i, :, k)) .* sign(mu_dks(:, k))' .* m_idks(i, :, k);

% Account for possible illegal values
if p_itk < 3 || p_itk > size(data, 1) - 2
    density = 0;
    return
end

if k == 1
    start = 1; 
else
    start = p_itks(i, t, k-1) + 3; % + 3
end

% Account for possible illegal values
if p_itk < start - 2 % p_itk < start
    density = 0;
    return
end

ends = p_itk - 3;
flat_sum = sum(data(start:ends, :).^2, 'all');

% Flat following bump k
start = p_itk + 3;

n = size(p_itks, 3);

if k == n
    ends = size(data, 1);
else
    ends = p_itks(i, t, k+1) - 3;
end

% Account for possible illegal values
if p_itk > ends + 2 % p_itk > ends
    density = 0;
    return
end

flat_sum = flat_sum + sum(data(start:ends, :).^2, 'all');

% Bump itself
bump_sum = sum( data(p_itk-2:p_itk+2, :).^2 - ... 
    2 * data(p_itk-2:p_itk+2, :) .* m_idks(i, :, k) .* P', 'all' );

% Experimentation with simple correlation
%bump_sum = -sum( data(p_itk-2:p_itk+2, :) .* m_idks(i, :, k) .* P', 'all' );
% bump_sum = -sum( data(p_itk-2:p_itk+2, :) .* (P' * m_idks(i, :, k)), 'all' );

density = exp(- (flat_sum + bump_sum) / V) * exp(-(p_itk^2 - 2*p_itk*nu_iks(i, k) + nu_iks(i, k)^2) / (2 * sigma2_iks(i, k)));
end

