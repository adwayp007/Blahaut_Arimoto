function [] = baa_dmc_cap_sim()
%% Calculates Binary Symmetric Channel (BSC) capacity using classical BAA

% Define the range of error probabilities
eps = 0:0.05:0.5;

% Initialize arrays to store channel capacity and entropy function values
channel_cap = zeros(1, length(eps));
entr_fun = zeros(1, length(eps));

% Iterate over each error probability
for i = 1:length(eps)
    disp(eps(i));
    
    % Initialize transition matrix
    tran_mat = [1 - eps(i), eps(i); eps(i), 1 - eps(i)];
    
    % Calculate channel capacity and input pmf using classical BAA
    [cap, in_pmf] = calculate_cap_dmc(tran_mat);
    
    % Store the calculated channel capacity and entropy function values
    channel_cap(i) = cap;
    entr_fun(i) = 1 - log2_entropy(eps(i), 1 / eps(i)) - log2_entropy(1 - eps(i), 1 / (1 - eps(i)));
end

% Plot the results
plot(eps, channel_cap, 'sr', eps, entr_fun, 'b', 'LineWidth', 1.5);
legend('Using BAA', '1-h(p)');
xlabel('Error Probability p');
ylabel('Capacity');
end
