function val = log2_entropy(A, B)
%% Calculate A .* log2(B) with handling of 0log0 scenario
%% Dimensions of A and B must match

% Calculate the element-wise product of A and the base-2 logarithm of B
val = A .* log2(B);

% Handle cases where log2(0) results in NaN by setting them to 0
val(isnan(val)) = 0;
end
