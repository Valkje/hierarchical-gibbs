function samples = invgamrnd(alpha, beta, n)
%INVGAMRND Summary of this function goes here
%   Detailed explanation goes here
arguments
   alpha double
   beta double
   n double = 1
end

    samples = 1 ./ gamrnd(alpha, 1 / beta, n, 1);
end

