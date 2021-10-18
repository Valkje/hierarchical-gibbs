function density = alpha_kPdf(alpha_k, c_k, d_k, e_k, sigma2_iks, beta_k, product)
%ALPHA_KPDF Summary of this function goes here
%   Detailed explanation goes here
    arguments
        alpha_k double
        c_k double
        d_k double
        e_k double
        sigma2_iks (:, 1) double
        beta_k double
        product double = 0
    end

    N = length(sigma2_iks);
    
    % Save some computational overhead if possible
    if product == 0
        density = beta_k^((e_k + N) * alpha_k) / ... 
            (gamma(alpha_k)^(d_k + N) * (c_k * prod(sigma2_iks))^(alpha_k + 1));
    else
        density = beta_k^((e_k + N) * alpha_k) / ... 
            (gamma(alpha_k)^(d_k + N) * (c_k * product)^(alpha_k + 1));
    end
end

