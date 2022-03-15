function visualize(params, hists)
%VISUALIZE Visualizes all parameters and their histories.
%   Creates bar plots, histograms and topoplots. Can become pretty slow for
%   many subjects with many bumps, so you might want to selectively comment
%   out some parts, or not call this function altogether when you need to
%   run many iterations.
    iter        = params.iter;
    x           = params.x;
    y           = params.y;
    subjects    = params.subjects;
    coeff10     = params.coeff10;
    latent10    = params.latent10;
    data        = params.data;
    chanlocs    = params.chanlocs;
    conds       = params.conds;
    cond        = params.cond;
    relative    = params.relative;
    firstVis    = params.firstVis;

    burnIn = params.burnIn;
    renderLag = params.renderLag;
    thinLag = params.thinLag;
    
    N = 1;
    n = size(hists.mu_dk, 3);
    D = 1;
    
    % Recorded samples to display
    displaySelection = burnIn+1:thinLag:iter;
    
    persistent h;
    if firstVis
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
    
    %% Topoplots
    set(0, 'CurrentFigure', h(10));
    clf(gcf)

    if isfield(params, 'n')
        n = params.n; % Incremental method, params.n is probably smaller than the corresponding histogram dimensions
    end
    
    % I use the jet colormap purely by convention, but you should know that
    % it has quite some deficiencies. Consider using the turbo colormap
    % instead.
    colors = jet(256);
    line_height = 0.05;
    
    nsubj = size(hists.m_idk, 2);

    trial_lens = y - x + 1;
    
    axes('Units','normalized','Position',[0 0 1 1]);
    ylim([0 1])
    xlim([0 1])
    
    % Plot every bump location as a line segment
    % If we don't do this here, we would have to redraw the axes for every
    % subject and lose the segments of the previous subject
    for i = 1:nsubj
        yline(1 - (i/nsubj))
        
        subject_lens = trial_lens(subjects == i & conds == cond);
        subject_trial_num = length(subject_lens);
        
        lines_x = squeeze(hists.p_itk(iter, i, 1:subject_trial_num, 1:n));
        
        % For relative bumps.
        if relative
            lines_x = cumsum(lines_x, 2);
        end
        
        lines_x = lines_x ./ subject_lens;
        
        ypos = 1 - (i/nsubj);
        
        
        for k = 1:n
            % Assume pc 1 is approx. in [-1.5, 1.5]. Bring up to positives
            pc = hists.m_idk(iter, i, 1, k) + 1.5;
            % Clamp pc 1 to [0, 3], then normalize to [0, 1]
            pc = min(3, max(0, pc)) / 3;
            % Use pc to get a color from colormap
            color = colors(max(1, min(ceil(pc * 256), 256)), :);
            
            lines_x_k = lines_x(:, k)';
            line([lines_x_k; lines_x_k], ... 
                [ypos; ypos+min(1/nsubj, line_height+k*0.02)], ... 
                'Color', color)
        end
    end

    for i = 1:nsubj
        subject_lens = trial_lens(subjects == i & conds == cond);
        subject_xpos = squeeze(hists.p_itk(iter, i, 1:length(subject_lens), 1:n));
        
        % For relative bumps
        if relative
            subject_xpos = cumsum(subject_xpos, 2);
        end
        
        subject_xpos = subject_xpos ./ subject_lens;
        subject_xpos = mean(subject_xpos, 1);
        
        ypos = 1 - (i/nsubj);
        
        % The width of the topologies is defined by the absolute bump width
        % (5 samples) over the mean trial length
        topo_size = 5 / mean(trial_lens(conds == cond));
        % On second thought, this is really small
        topo_size = topo_size * 2;
        topo_size = min(1/nsubj, topo_size);
        
        % Plot the topologies
        topo = reconstruct(squeeze(hists.m_idk(iter, i, :, 1:n)), coeff10, latent10, mean(data));
        for k = 1:n
            % xpos should correspond to centre of topo image
            xpos = subject_xpos(k) - topo_size / 2;
            
            axes('Units','normalized','Position',[xpos ypos topo_size topo_size]);
            topoplot(topo(k,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');
        end
    end
    
    set(0, 'CurrentFigure', h(11));
    clf(gcf)
    for k = 1:n
        subplot(2, round(n/2), k)
        histogram(hists.p_itk(displaySelection, 1, 1, k))
        title(['p_{itk}, i=1, t=1, k=' num2str(k)])
        set(gca, 'FontSize', 14)
    end
    
    drawnow
end

