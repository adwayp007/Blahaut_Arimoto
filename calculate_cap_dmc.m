function [cap,in_pmf] = calculate_cap_dmc(W)

[m,n]=size(W);

q_0 = rand(1,m);
q_0 = q_0/sum(q_0);

q_k = q_0;
eps = 0.0001;
error=1;
for i=1:1000
    t_fun = zeros(1,m);
    r = q_k*W;
    for j=1:m
        for k=1:n
            temp = W(j,k)*log2(q_k(1,j)*W(j,k)/r(1,k));
            t_fun(1,j)=t_fun(1,j)+temp;
        end
        
    end
    q_k = 2.^(t_fun);
    sum1 = sum(q_k);
    q_k = q_k/sum1;
    error = max(t_fun - log2(q_k))-min(t_fun - log2(q_k));
    
end

in_pmf = q_k;

cap = sum(q_k.*(log2(1./q_k)+t_fun));

end


