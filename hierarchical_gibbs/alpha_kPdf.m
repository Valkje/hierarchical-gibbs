function density = alpha_kPdf(alpha_k, a_k, b_k, q_k, alpha_iks, beta_k, product)
%ALPHA_KPDF Summary of this function goes here
%   Detailed explanation goes here
    arguments
        alpha_k double
        a_k double
        b_k double
        q_k double
        alpha_iks (:, 1) double
        beta_k double
        product double = 0
    end

    N = length(alpha_iks);
    
    if product == 0
        density = beta_k^((q_k + N) * alpha_k) / ... 
            (gamma(alpha_k)^(b_k + N) * (a_k * prod(alpha_iks))^(alpha_k + 1));
    else
        density = beta_k^((q_k + N) * alpha_k) / ... 
            (gamma(alpha_k)^(b_k + N) * (a_k * product)^(alpha_k + 1));
    end
end

