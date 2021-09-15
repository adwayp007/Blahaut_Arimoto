function [val] = log2_entropy(A,B)
%% Calculates A.*log_b(B) taking care of 0log 0 scenario
%% dimensions of A nd B must match

val = A.*log2(B);

val(isnan(val))=0;

end

