%% Test EC LQ

function [dX,dU,dL,dM,dP,sn] = solve_ec_single_game(F,... % time-indexed dict of dynamics
                                               H,... % time-indexed dict of player-indexed constraints
                                               Q,... % time-indexed dict of player-indexed costs
                                               N,... % number of players
                                               T,... % number of timesteps
                                               KK,... % precomputed policies
                                               kk,... % precomputed policies
                                               p,... % ego player index
                                               x0)   % initial state
    
    
    n = size(F{1},1);
    all_m = 0;
    for i = 1:N
        start_m{i} = all_m;
        m{i} = size(H{1,i},2) - 1 - n;
        all_m = all_m + m{i};
    end
    all_l = 0;
    for t = 1:T
        l{t} = size(H{t,p},1);
        H{t,p} = [H{t,p}(:,1:1+n) zeros(l{t},start_m{p}) H{t,p}(:,2+n:end) zeros(l{t},all_m-m{p})];
        all_l = all_l+l{t};
    end
    l{T+1} = size(H{T+1,p},1);
    all_l = all_l+l{T+1};
    
    num_primals = (T+1)*n + T*all_m;
    num_dynamic_constraints = (T+1)*n;
    num_other_constraints = all_l;
    num_policy_constraints = T*(all_m-m{p});
    
    all_vars = num_primals+num_dynamic_constraints+num_other_constraints + num_policy_constraints;
    
    M = zeros(all_vars);
    b = zeros(all_vars,1);
    
    ind = 0;
    for t = 1:T
        M(ind+1:ind+n+all_m,ind+1:ind+n+all_m) = 0.5*Q{t,p}(2:end,2:end); % trick used for easy filling 
        b(ind+1:ind+n+all_m) = -Q{t,p}(2:end,1);
        ind = ind+n+all_m;
    end
    M(ind+1:ind+n,ind+1:ind+n) = 0.5*Q{T+1,p}(2:end,2:end);
    b(ind+1:ind+n) = -Q{T+1,p}(2:end,1);
    ind = ind+n;
    
    ind2 = 0;
    M(ind+1:ind+n,1:n) = -eye(n);
    b(ind+1:ind+n) = x0;
    ind = ind+n;
    for t = 1:T
        b(ind+1:ind+n) = -F{t}(:,1);
        M(ind+1:ind+n,ind2+1:ind2+n+all_m) = F{t}(:,2:end);
        M(ind+1:ind+n,ind2+n+all_m+1:ind2+n+all_m+n) = -eye(n);
        ind2 = ind2+n+all_m;
        ind = ind+n;
    end
   
    ind2 = 0;
    for t = 1:T
        M(ind+1:ind+l{t},ind2+1:ind2+n+all_m) = H{t,p}(:,2:end);
        b(ind+1:ind+l{t}) = -H{t,p}(:,1);
        ind = ind+l{t};
        ind2 = ind2+n+all_m;
    end
    M(ind+1:ind+l{T+1},ind2+1:ind2+n) = H{T+1,p}(:,2:end);
    b(ind+1:ind+l{T+1}) = -H{T+1,p}(:,1);
    ind = ind+l{T+1};
    
    ind2 = 0;
    for t = 1:T
        fI = eye(all_m-m{p});
        M(ind+1:ind+all_m-m{p},ind2+1:ind2+n) = KK{t,p};
        M(ind+1:ind+all_m-m{p},ind2+n+1:ind2+n+start_m{p}) = -fI(:,1:start_m{p});
        M(ind+1:ind+all_m-m{p},ind2+n+start_m{p}+m{p}+1:ind2+n+all_m) = -fI(:,start_m{p}+1:end);
        b(ind+1:ind+all_m-m{p}) = -kk{t,p};
        ind = ind+all_m-m{p};
        ind2 = ind2 + n+all_m;
    end
    M = M+M';
    sol = M\b;
    sn = norm(sol(1:num_primals));
    
    dX = cell(T+1,1);
    dU = cell(T,N);
    dL = cell(T+1,N);
    dM = cell(T+1,N);
    dP = cell(T,N);
    
    ind = 0;
    for t = 1:T
        dX{t} = sol(ind+1:ind+n);
        ind = ind+n;
        for i = 1:N
           dU{t,i} = sol(ind+1:ind+m{i});
           ind = ind+m{i};
        end
    end
    dX{T+1} = sol(ind+1:ind+n);
    ind = ind+n;
    for t = 1:T+1
        dL{t,p} = sol(ind+1:ind+n);
        ind = ind+n;
    end
    for t = 1:T+1
        dM{t,p} = sol(ind+1:ind+l{t});
        ind = ind+l{t};
    end
    for t = 1:T
        dP{t,p} = sol(ind+1:ind+all_m-m{p});
        ind = ind+all_m-m{p};
    end
end