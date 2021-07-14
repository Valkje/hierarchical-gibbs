function visualize(max_iter, iter, ... 
    normedscore10, x, y, subjects, ... 
    coeff10, latent10, data, chanlocs, ...
    alpha_k_hist, beta_k_hist, alpha_ik_hist, t_i_hist, ... 
    mu_dk_hist, tau2_dk_hist, m_idk_hist)
%VISUALIZE Visualizes all parameters and their histories.
%   Creates bar plots, histograms and topoplots.
    N = 1;
    n = size(mu_dk_hist, 3);
    D = 1;
    
    burnIn = 1000;
    renderLag = 1000;
    thinLag = 20;
    
    % Recorded samples to display
    displaySelection = burnIn+1:thinLag:iter;
    
    persistent h;
    if iter == 1
        numFigs = 10;
        h = gobjects(numFigs, 1);
        for i = 1:numFigs
            h(i) = figure(i);
            clf(h(i))
        end
    end

    if mod(iter, renderLag) ~= 0 && iter ~= 1
        return
    end
    
    set(0, 'CurrentFigure', h(1));
    clf(gcf)
    barh(squeeze(t_i_hist(iter, :, :)), 'stacked')
    legs = cell(n+1, 1);
    for k = 1:n+1
        legs{k} = ['Flat ' num2str(k)];
    end
    set(gca, 'FontSize', 14)
    
    title('Stacked duration proportions per participant')
    xlabel('t_i')
    ylabel('Participant')
    
    legend(legs)
    
    set(0, 'CurrentFigure', h(2));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(t_i_hist(displaySelection, N, k))
        title(['t_{ik}, i=' num2str(N) ', k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(3));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(alpha_ik_hist(displaySelection, N, k))
        title(['\alpha_{ik}, i=' num2str(N) ', k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(4));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(alpha_k_hist(displaySelection, k))
        title(['\alpha_{k}, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(5));
    clf(gcf)
    for k = 1:n+1
        subplot(2, round((n+1)/2), k)
        histogram(beta_k_hist(displaySelection, k))
        title(['\beta_{k}, k=' num2str(k)])
    end
    
    set(0, 'CurrentFigure', h(6));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(m_idk_hist(displaySelection, N, d, k))
            title(['m_{idk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(7));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(mu_dk_hist(displaySelection, d, k))
            title(['\mu_{dk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(8));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(tau2_dk_hist(displaySelection, d, k))
            title(['\tau^2_{dk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(9));
    clf(gcf)
    
    nsubj = 10;
    
    bump_positions = cumsum(squeeze(t_i_hist(iter, :, :)), 2);
    
    for i = 1:nsubj 
        
        topo = reconstruct(squeeze(m_idk_hist(iter, i, :, :)), coeff10, latent10, mean(data));
        for k = 1:n
            ypos = 1 - (i/nsubj);
            xpos = bump_positions(i, k) / 1.1;

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
                    [starts, ends] = getFlatIndices(x, y, n, t_i_hist(iter, 1, :));

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
end

