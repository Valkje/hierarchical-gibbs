function visualize(max_iter, iter, ... 
    normedscore10, x, y, subjects, ... 
    coeff10, latent10, data, chanlocs, ...
    alpha_ks, beta_ks, alpha_iks, t_is, mu_dks, tau2_dks, m_idks)
%VISUALIZE Summary of this function goes here
%   Detailed explanation goes here
    N = 1;
    n = size(mu_dks, 2);
    D = 1;
    
    burnIn = 3000;
    lag = 100;
    
    persistent alpha_k_hist;
    persistent beta_k_hist;
    persistent alpha_ik_hist;
    persistent t_is_hist;
    persistent m_idks_hist;
    persistent mu_dks_hist;
    persistent tau2_dks_hist;
    persistent h;
    if iter == 1
        alpha_k_hist = zeros(max_iter, n+1);
        beta_k_hist = zeros(max_iter, n+1);
        alpha_ik_hist = zeros(max_iter, n+1);
        t_is_hist = zeros(max_iter, n+1);
        m_idks_hist = zeros(max_iter, D, n);
        mu_dks_hist = zeros(max_iter, D, n);
        tau2_dks_hist = zeros(max_iter, D, n);
        
        numFigs = 10;
        h = gobjects(numFigs, 1);
        for i = 1:numFigs
            h(i) = figure(i);
        end
    end
    
    iter = iter + 1;
    
    if iter-burnIn > 0
        alpha_k_hist(iter-burnIn, :) = alpha_ks;
        beta_k_hist(iter-burnIn, :) = beta_ks;
        alpha_ik_hist(iter-burnIn, :) = alpha_iks(N, :);
        t_is_hist(iter-burnIn, :) = t_is(N, :);
        m_idks_hist(iter-burnIn, :, :) = m_idks(N, 1:D, :);
        mu_dks_hist(iter-burnIn, :, :) = mu_dks(1:D, :);
        tau2_dks_hist(iter-burnIn, :, :) = tau2_dks(1:D, :);
    end

    if mod(iter, lag) ~= 0 && iter ~= 1
        return
    end
    
    set(0, 'CurrentFigure', h(1));
    clf(gcf)
    barh(t_is, 'stacked')
    legs = cell(n+1, 1);
    for k = 1:n+1
        legs{k} = ['Flat ' num2str(k)];
    end
    
    title('Stacked duration proportions per participant')
    xlabel('t_i')
    ylabel('Participant')
    
    legend(legs)
    
    set(0, 'CurrentFigure', h(2));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(t_is_hist(1:iter-burnIn, k))
        title(['t_{ik}, i=' num2str(N) ', k=' num2str(k)])
    end
    
    set(0, 'CurrentFigure', h(3));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(alpha_ik_hist(1:iter-burnIn, k))
        title(['\alpha_{ik}, i=' num2str(N) ', k=' num2str(k)])
    end
    
    set(0, 'CurrentFigure', h(4));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(alpha_k_hist(1:iter-burnIn, k))
        title(['\alpha_{k}, k=' num2str(k)])
    end
    
    set(0, 'CurrentFigure', h(5));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(beta_k_hist(1:iter-burnIn, k))
        title(['\beta_{k}, k=' num2str(k)])
    end
    
    set(0, 'CurrentFigure', h(6));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(m_idks_hist(1:iter-burnIn, d, k))
            title(['m_{idk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
    end
    
    set(0, 'CurrentFigure', h(7));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(mu_dks_hist(1:iter-burnIn, d, k))
            title(['\mu_{dk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
    end
    
    set(0, 'CurrentFigure', h(8));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(tau2_dks_hist(1:iter-burnIn, d, k))
            title(['\tau^2_{dk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
    end
    
    set(0, 'CurrentFigure', h(9));
    clf(gcf)
    
    nsubj = 2;
    y_offset = 0.5 / nsubj;
    
    for i = 1:nsubj 
        topo = reconstruct(squeeze(m_idks(i, :, :)), coeff10, latent10, mean(data));
        for k = 1:n
            xpos = (k - 1) / n;
            ypos = 1 - i / nsubj + y_offset;

            axes('Units','normalized','Position',[xpos ypos .1 .1]);
            topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');
        end
    end
    
    set(0, 'CurrentFigure', h(10));
    clf(gcf)
    
    if n == 2
        x = x(subjects == 1);
        y = y(subjects == 1);
        
        densities = zeros(101, 101);
        i = 0;

        for t1 = 0:0.01:1
            i = i + 1;
            j = 0;

            for t2 = 0:0.01:1-t1
                j = j + 1;

                t3 = 1 - t1 - t2;

                ts = [t1 t2 t3];

                product = 1;
                for k = 1:3
                    [starts, ends] = getFlatIndices(x, y, n, t_is(1, :));

                    S_sqr_sum = 0;
                    for l = 1:length(x)
                        S_sqr_sum = S_sqr_sum + sum( normedscore10(starts(l):ends(l)).^2 );
                    end
                    
                    dens = log(ts(k) ^ (alpha_iks(1, k) - 1)) + -S_sqr_sum / 5;
                    
                    if ~isreal(dens)
                        dens = 0;
                    end
                    
                    if dens ~= -Inf
                        product = product + dens;
                    end
                end

                densities(i, j) = product;
            end
        end
        
        densities(densities ~= 0) = exp(densities(densities ~= 0) - median(densities(densities ~= 0)));
%         densities(densities > 5 & densities ~= Inf) = 5;

        surf(densities')
    end
    
    drawnow
    
%     figure(2)
%     clf(gcf)
%     hold on
%     
%     subj = 1;
%     ts = 0:0.01:1;
%     for k = 1:size(alpha_iks, 2)
%         plot(ts, ts.^(alpha_iks(subj, k)) - 1);
%     end
%     hold off

    % Distribution for alpha_iks
%     figure(3)
%     clf(gcf)
%     hold on
%     
%     x = 0:0.1:20;
%     for k = 1:length(alpha_ks)
%         densities = beta_ks(k) ^ alpha_ks(k) .* x .^ (-alpha_ks(k) - 1) .* exp(-beta_ks(k) ./ x) ./ gamma(alpha_ks(k));
%         plot(x, densities, 'DisplayName', num2str(k))
%     end
%     legend
%     
%     hold off
end

