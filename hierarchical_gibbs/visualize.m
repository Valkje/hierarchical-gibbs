function visualize(max_iter, iter, start_iter,... 
    normedscore10, x, y, subjects, ... 
    coeff10, latent10, data, chanlocs, ...
    hists)
%VISUALIZE Visualizes all parameters and their histories.
%   Creates bar plots, histograms and topoplots.
    N = 1;
    n = size(hists.mu_dk, 3);
    D = 1;
    
    burnIn = 10;
    renderLag = 10;
    thinLag = 1;
    
    % Recorded samples to display
    displaySelection = burnIn+1:thinLag:iter;
    
    persistent h;
    if iter == start_iter
        numFigs = 11;
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
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.nu_ik(displaySelection, N, k))
        title(['\nu_{ik}, i=' num2str(N) ', k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(2));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.sigma2_ik(displaySelection, N, k))
        title(['\sigma^2_{ik}, i=' num2str(N) ', k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(3));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.xi_k(displaySelection, k))
        title(['\xi_{k}, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(4));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.upsilon2_k(displaySelection, k))
        title(['\upsilon^2_{k}, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(5));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.alpha_k(displaySelection, k))
        title(['\alpha_{k}, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(6));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.beta_k(displaySelection, k))
        title(['\beta_{k}, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
%     set(0, 'CurrentFigure', h(1));
%     clf(gcf)
%     barh(squeeze(t_i_hist(iter, :, :)), 'stacked')
%     legs = cell(n+1, 1);
%     for k = 1:n+1
%         legs{k} = ['Flat ' num2str(k)];
%     end
%     set(gca, 'FontSize', 14)
%     
%     title('Stacked duration proportions per participant')
%     xlabel('t_i')
%     ylabel('Participant')
%     
%     legend(legs)
    
%     set(0, 'CurrentFigure', h(2));
%     clf(gcf)
%     for k = 1:n+1
%         subplot(2, round((n+1)/2), k)
%         histogram(t_i_hist(displaySelection, N, k))
%         title(['t_{ik}, i=' num2str(N) ', k=' num2str(k)])
%         set(gca, 'FontSize', 14)
%     end
%     
%     set(0, 'CurrentFigure', h(3));
%     clf(gcf)
%     for k = 1:n+1
%         subplot(2, round((n+1)/2), k)
%         histogram(alpha_ik_hist(displaySelection, N, k))
%         title(['\alpha_{ik}, i=' num2str(N) ', k=' num2str(k)])
%         set(gca, 'FontSize', 14)
%     end
%     
%     set(0, 'CurrentFigure', h(4));
%     clf(gcf)
%     for k = 1:n+1
%         subplot(2, round((n+1)/2), k)
%         histogram(alpha_k_hist(displaySelection, k))
%         title(['\alpha_{k}, k=' num2str(k)])
%         set(gca, 'FontSize', 14)
%     end
%     
%     set(0, 'CurrentFigure', h(5));
%     clf(gcf)
%     for k = 1:n+1
%         subplot(2, round((n+1)/2), k)
%         histogram(beta_k_hist(displaySelection, k))
%         title(['\beta_{k}, k=' num2str(k)])
%     end
    
    set(0, 'CurrentFigure', h(7));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(hists.m_idk(displaySelection, N, d, k))
            title(['m_{idk}, d=' num2str(d) ', k=' num2str(k)])
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
            histogram(hists.mu_dk(displaySelection, d, k))
            title(['\mu_{dk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(9));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        hold on
        for d = 1:D
            histogram(hists.tau2_dk(displaySelection, d, k))
            title(['\tau^2_{dk}, d=' num2str(d) ', k=' num2str(k)])
        end
        hold off
        set(gca, 'FontSize', 14)
    end
    
    set(0, 'CurrentFigure', h(10));
    clf(gcf)

    nsubj = 2;

%         bump_positions = cumsum(squeeze(t_i_hist(iter, :, :)), 2);

    for i = 1:nsubj 

        topo = reconstruct(squeeze(hists.m_idk(iter, i, :, :)), coeff10, latent10, mean(data));
        for k = 1:n
            ypos = 1 - (i/nsubj);
            xpos = k * (1 / (n + 1));

            axes('Units','normalized','Position',[xpos ypos .1 .1]);
            topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');
        end
    end
    
    set(0, 'CurrentFigure', h(11));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.p(displaySelection, k))
        title(['p_{itk}, i=1, t=1, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    drawnow
end

