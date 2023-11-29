function [] = bsc_capacity_rll_source()
%% Source cannot output more than one consecutive 0
%% Channel is Binary Symmetric Channel (BSC)

% Generate random probabilities for the transition matrix initialization
r = rand(1, 2);
q = [0, r(1, 1); r(1, 1), r(1, 2)] / (sum(r) + r(1, 1)); % Initialization

% Define error probabilities
eps = 0:0.05:0.5;

% Stopping criteria threshold
delta = 0.001;

% Initialize result arrays
q_r = q;
cap_results = zeros(1, length(eps));
entr_fun = zeros(1, length(eps));

% Iterate over each error probability
for iter = 1:length(eps)
    disp(eps(1, iter));
    
    % Initialize capacity variables
    cap = 0;
    prev_cap = 1;

    % Iterate until convergence or cap<0
    while (abs(cap - prev_cap) > delta || cap < 0)
        % Calculate parameters for the transition matrix
        T = calculate_T_params_dmc(q_r, eps(1, iter));
        T(isnan(T)) = 0;
        A = 2.^(T);
        A(1, 1) = 0;

        % Eigenvalue decomposition
        [M, L] = eig(A);
        [m, index] = max(diag(L));
        eig_vec = M(:, index);

        % Update transition matrix P
        P = zeros(2, 2);
        for i = 1:2
            for j = 1:2
                P(i, j) = (eig_vec(j, 1) / eig_vec(i, 1)) * (A(i, j) / m);
            end
        end

        % Update mu
        mu = [P(2, 1) / (P(2, 1) + P(1, 2)), P(1, 2) / (P(2, 1) + P(1, 2))];

        % Update q_r using P and mu
        for i = 1:2
            for j = 1:2
                q_r(i, j) = P(i, j) * mu(1, i);
            end
        end

        % Update capacity
        prev_cap = cap;
        cap = sum(sum(log2_entropy(q_r, 1./P) + q_r.*T));
        
        disp("The capacity at this stage is : " + cap);
        cap_results(1, iter) = cap;
    end

    % Calculate entropy function
    entr_fun(1, iter) = 1 - log2_entropy(eps(1, iter), 1 / eps(1, iter)) - log2_entropy(1 - eps(1, iter), 1 / (1 - eps(1, iter)));
end

% Plot the results
plot(eps, cap_results, 'b', 'LineWidth', 1.4);
hold on;
plot(eps, entr_fun, 'g', 'LineWidth', 1.4);
legend('Capacity of BSC using RLL source', 'Capacity of BSC in general');
xlabel('Error prob');
ylabel('Capacity');

% End of the function
end
