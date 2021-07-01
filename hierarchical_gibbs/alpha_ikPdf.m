function density = alpha_ikPdf(alpha_ik, alpha_iks, t_ik, alpha_k, beta_k)
%ALPHA_IKPDF Summary of this function goes here
%   Detailed explanation goes here
    density = (gamma(sum(alpha_iks) + alpha_ik) * t_ik^(alpha_ik - 1)) / ...
        (gamma(alpha_ik) * alpha_ik^(alpha_k + 1)) * ... 
        exp(-beta_k / alpha_ik);
end

