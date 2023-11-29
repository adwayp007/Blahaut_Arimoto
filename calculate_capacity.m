function [I] = calculate_capacity(A)
%% Calculate channel capacity

% Get the size of matrix A
[n, m] = size(A);

% Initialize the channel capacity
I = 0;

% Iterate over each column of matrix A
for i = 1:m
    % Calculate the conditional probability distribution
    p_x_given_y = A(:, i) / sum(A(:, i) * 0.5);

    % Handle cases where log2(0) results in -Inf
    temp = log2(p_x_given_y);
    temp(isinf(temp)) = 0;

    % Update the channel capacity
    I = I + 0.5 * sum(A(:, i) .* temp);
end
end
