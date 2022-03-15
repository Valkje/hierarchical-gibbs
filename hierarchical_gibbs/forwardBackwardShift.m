function [p_itks] = forwardBackwardShift(p_itks, subjects, conds, ...
        trial_lens, trial_nums, N, n, cond)
%FORWARDBACKWARDSHIFT Shifts bumps such that they do not (fully) overlap.
%   Some bump initialization methods might place two bumps in the same
%   trial at the same position. Although partial overlap is allowed,
%   complete overlap is not. This function first shifts all overlapping
%   bumps towards the end of the trial. Once that has been done, no bumps
%   overlap anymore, but it is possible that some extend beyond the trial
%   border. To remedy that problem, such bumps are shifted backwards; if
%   that results in new overlap, then the offending bumps are shifted
%   backwards as well.
    
    for i = 1:N
        subject_lens = trial_lens(subjects == i & conds == cond);

        for t = 1:trial_nums(i)
            trial_len = subject_lens(t);

            %% Forward shift
            
            % A bump is five indices wide, so it cannot be placed at a
            % position smaller than 3 (because then it would extend beyond
            % the start of the trial).
            if p_itks(i, t, 1) < 3
                p_itks(i, t, 1) = 3;
            end

            % For every bump k that is at or below the location of bump
            % k-1, place it just after bump k-1.
            for k = 2:n
                if p_itks(i, t, k) < p_itks(i, t, k-1) + 1
                    p_itks(i, t, k) = p_itks(i, t, k-1) + 1;
                end
            end

            %% Backward shift
            
            % A bump is five indices wide, so it cannot be placed at a
            % position larger than trial_len-3 (because then it would
            % extend beyond the end of the trial).
            if p_itks(i, t, end) > trial_len - 3
                p_itks(i, t, end) = trial_len - 3;
            end

            % For every bump k that is at or above the location of bump
            % k+1, place it just before bump k+1.
            for k = n-1:-1:1
                if p_itks(i, t, k) > p_itks(i, t, k+1) - 1
                    p_itks(i, t, k) = p_itks(i, t, k+1) - 1;
                end
            end
        end
    end
    
end

