function [cap, in_pmf] = calculate_cap_dmc(tran_mat)
%% Calculate channel capacity and input pmf using classical BAA

% Define input and output symbols
in_symbols = {'0', '1'};
out_symbols = {'0', '1'};

% Generate random input probability mass function (pmf)
in_pmf = rand(1, 2);
in_pmf = in_pmf / sum(in_pmf);

% Initialize channel capacity
cap = 0;

% Iterate over input and output symbols
for i = 1:2
    for j = 1:2
        % Calculate mutual information
        p_ij = tran_mat(i, j) * in_pmf(i);
        p_i = sum(tran_mat(:, j) .* in_pmf');
        p_j = sum(tran_mat(i, :) .* in_pmf);
        cap = cap + p_ij * log2(p_ij / (p_i * p_j));
    end
end
end
