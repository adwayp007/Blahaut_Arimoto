function [] = baa_dmc_cap_sim()

%% Calculates BSC capacity using classical BAA
eps=0:0.05:0.5;

channel_cap = zeros(1,length(eps));
entr_fun = zeros(1,length(eps));
for i=1:length(eps)
    
    tran_mat = [1-eps(1,i) eps(1,i);eps(1,i) 1-eps(1,i)];
    [cap, in_pmf] = calculate_cap_dmc(tran_mat);
    channel_cap(1,i) = cap;
    entr_fun(1,i) = 1-log2_entropy(eps(1,i),1/eps(1,i))-log2_entropy(1-eps(1,i),1/(1-eps(1,i)));
    
end

plot(eps,channel_cap,'sr',eps,entr_fun,'b','LineWidth',1.5);
legend('Using BAA','1-h(p)');
xlabel('Error Probability p');
ylabel('Capacity');