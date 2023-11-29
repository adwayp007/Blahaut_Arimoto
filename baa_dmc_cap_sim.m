function [] = baa_dmc_cap_sim()
%% Calculates BSC capacity using classical BAA

% Define the range of error probabilities
eps = 0:0.05:0.5;

% Initialize arrays to store channel capacity and entropy function values
channel_cap = zeros(1, length(eps));
entr_fun = zeros(1, length(eps));

% Iterate over each error probability
for i = 1:length(eps)
    
    % Transition matrix for the Binary Symmetric Channel (BSC)
    tran_mat = [1 - eps(1, i), eps(1, i); eps(1, i), 1 - eps(1, i)];
    
    % Calculate channel capacity and input pmf using classical BAA
    [cap, in_pmf] = calculate_cap_dmc(tran_mat);
    
    % Store the calculated channel capacity and entropy function values
    channel_cap(1, i) = cap;
    entr_fun(1, i) = 1 - log2_entropy(eps(1, i), 1 / eps(1, i)) - log2_entropy(1 - eps(1, i), 1 / (1 - eps(1, i)));
    
end

% Plot the results
plot(eps, channel_cap, 'sr', eps, entr_fun, 'b', 'LineWidth', 1.5);
legend('Using BAA', '1-h(p)');
xlabel('Error Probability p');
ylabel('Capacity');

% End of the function
end
