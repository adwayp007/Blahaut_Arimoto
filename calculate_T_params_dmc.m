function [T_hat] = calculate_T_params_dmc(Q,eps)

    N=100000;
    s_states = {'0','1'};
    y_states = {'0','1'};

    mu = sum(Q);
    T=zeros(2,2);
    for i=1:2
        for j=1:2
            T(i,j) = Q(i,j)/mu(1,i);  %% Source state transition matrix
        end
    end

    Emis = 0.5*[(1-eps)+(1-eps) ...
        (eps)+(eps); ...
        (eps)+(eps) ...
        (1-eps)+(1-eps)];

    [out_seq, int_s_states] = hmmgenerate(N+1,T,Emis,'Symbols',...
        y_states,'Statenames',s_states);

    b = cell(1,N);
    for i=1:N
        temp1 = int_s_states{i};
        temp2 = int_s_states{i+1};
        b(1,i) = cellstr(strcat(temp1,temp2));
    end

    V_s_y = hmmdecode(out_seq(2:end),T,Emis,'Symbols',{'0','1'});

    T_branch = [T(1,:) 0 0; 0 0 T(2,:);T(1,:) 0 0; 0 0 T(2,:)];

    Emis_branch = [Emis;Emis];

    V_b_y = hmmdecode(out_seq(2:end),T_branch,Emis_branch,...
        'Symbols',{'0','1'});

    T_hat = zeros(2,2);

    T_hat(1,1) = (sum(log2_entropy(V_b_y(1,:),V_b_y(1,:)))/Q(1,1) ...
        - sum(log2_entropy(V_s_y(1,:),V_s_y(1,:)))/mu(1,1))/N;
    T_hat(1,2) = (sum(log2_entropy(V_b_y(2,:),V_b_y(2,:)))/Q(1,2) ...
        - sum(log2_entropy(V_s_y(1,:),V_s_y(1,:)))/mu(1,1))/N;
    T_hat(2,1) = (sum(log2_entropy(V_b_y(3,:),V_b_y(3,:)))/Q(2,1) ...
        - sum(log2_entropy(V_s_y(2,:),V_s_y(2,:)))/mu(1,2))/N;
    T_hat(2,2) = (sum(log2_entropy(V_b_y(4,:),V_b_y(4,:)))/Q(2,2) ...
        - sum(log2_entropy(V_s_y(2,:),V_s_y(2,:)))/mu(1,2))/N;


end

