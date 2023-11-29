function [T_hat] = calculate_T_params_dmc(Q, eps)
%% Calculate parameters for transition matrix T using Hidden Markov Model (HMM)

% Set the number of samples
N = 100000;

% Define state and symbol sets
s_states = {'0', '1'};
y_states = {'0', '1'};

% Compute mu from the given transition matrix Q
mu = sum(Q);

% Build the source state transition matrix T
T = Q ./ mu;

% Define emission matrix Emis for Binary Symmetric Channel (BSC)
Emis = 0.5 * [(1 - eps) + (1 - eps), (eps) + (eps); (eps) + (eps), (1 - eps) + (1 - eps)];

% Generate an output sequence using HMM
[out_seq, int_s_states] = hmmgenerate(N + 1, T, Emis, 'Symbols', y_states, 'Statenames', s_states);

% Create state pairs
b = cell(1, N);
for i = 1:N
    b{i} = strcat(int_s_states{i}, int_s_states{i + 1});
end

% Decode sequences using HMM
V_s_y = hmmdecode(out_seq(2:end), T, Emis, 'Symbols', {'0', '1'});

% Expand transition matrix and emission matrix for branching HMM
T_branch = [T(1, :) 0 0; 0 0 T(2, :); T(1, :) 0 0; 0 0 T(2, :)];

Emis_branch = repmat(Emis, 2, 1);

V_b_y = hmmdecode(out_seq(2:end), T_branch, Emis_branch, 'Symbols', {'0', '1'});

% Calculate the estimated transition matrix T_hat
T_hat = zeros(2, 2);
for i = 1:2
    for j = 1:2
        idx = 2 * (i - 1) + j;
        T_hat(i, j) = (sum(log2_entropy(V_b_y(idx, :), V_b_y(idx, :))) / Q(i, j) ...
            - sum(log2_entropy(V_s_y(i, :), V_s_y(i, :))) / mu(i)) / N;
    end
end
end
