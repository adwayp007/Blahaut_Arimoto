function [I] = calculate_capacity(A)

[n,m] = size(A);
I =0;
for i=1:m
    temp=log2(A(:,i)./sum(A(:,i)*0.5));
    temp(isinf(temp))=0;
    I = I+0.5*sum((A(:,i).*(temp)));
end
end

