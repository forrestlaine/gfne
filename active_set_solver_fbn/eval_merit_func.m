function merit = eval_merit_func(terms,K,N,T,m,start_m)
    x = terms{1};
    u = terms{2};
    lambda = terms{3};
    mu = terms{4};
    gam = terms{5};
    psi = terms{6};
    F = terms{7};
    H = terms{8};
    G = terms{9};
    Q = terms{10};
    all_m = 0;
    for i = 1:N
        all_m = all_m + m{i};
    end

    n = numel(x{1});
    merit = 0;
    
    for i = 1:N
        for t = 1:T
            primal  = [1;x{t}];
            
            for j = 1:N
                primal = [primal; u{t,j}];
            end

                if t == 1
                    lag = Q{t,i}(n+2+start_m{i}:n+1+start_m{i}+m{i},1) + ...
                          F{t}(:,n+2+start_m{i}:n+1+start_m{i}+m{i})'*lambda{t,i} - ...
                          H{t,i}(:,n+2:end)'*mu{t,i} - ...
                          G{t,i}(:,n+2:end)'*gam{t,i};
                else
                    KK = [];
                    HH = [H{t,i}(:,1:n+1)];
                    GG = [G{t,i}(:,1:n+1)];
                    for j = 1:N
                        if i~=j
                            KK = [KK; K{t,j}];
                            HH = [HH zeros(size(HH,1),m{j})];
                            GG = [GG zeros(size(GG,1),m{j})];
                        else
                            HH = [HH H{t,i}(:,n+2:end)];
                            GG = [GG G{t,i}(:,n+2:end)];
                        end
                    end
                    EI = eye(all_m);
                    EI = [EI(1:start_m{i},:); EI(start_m{i}+m{i}+1:end,:)];
                    KK = [KK -EI];
                    lag = Q{t,i}(2:end,1) + ...
                          F{t}(:,2:end)'*lambda{t,i} - ...
                          [lambda{t-1,i}; zeros(all_m,1)] - ...
                          HH(:,2:end)'*mu{t,i} - ...
                          GG(:,2:end)'*gam{t,i} + ...
                          KK'*psi{t-1,i};
                end
                resid_eq = H{t,i}(:,1);
                resid_ineq = min(G{t,i}(:,1),0);
                resid_ineq_mult = min(gam{t,i},0);
                resid_dyn = F{t}(:,1);
                resid_comp = G{t,i}(:,1)'*gam{t,i};
                cur_merit = norm(lag) + ...
                                norm(resid_eq) + ...
                                norm(resid_ineq) + ...
                                norm(resid_ineq_mult) + ...
                                norm(resid_dyn) + ...
                                norm(resid_comp);
                merit = merit + cur_merit;
        end
           
        primal = [1;x{T+1}];
        lag = Q{T+1,i}(2:end,1) - ...
              lambda{T,i}- ...
              H{T+1,i}(:,2:end)'*mu{T+1,i} - ...
              G{T+1,i}(:,2:end)'*gam{T+1,i};
        resid_eq = H{T+1,i}(:,1);
        resid_ineq = min(G{T+1,i}(:,1),0);
        resid_ineq_mult = min(gam{T+1,i},0);
        resid_comp = G{T+1,i}(:,1)'*gam{T+1,i};
        merit = merit + norm(lag) + ...
                        norm(resid_eq) + ...
                        norm(resid_ineq) + ...
                        norm(resid_ineq_mult) + ...
                        norm(resid_comp);
    end
      
end