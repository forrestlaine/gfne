%% Test EC LQ

function [K,k,L,l] = solve_ec_lq_game_b(F,... % time-indexed dict of dynamics
                                        H,... % time-indexed dict of player-indexed constraints
                                        Q,... % time-indexed dict of player-indexed costs
                                        N,... % number of players
                                        T)     
    
    K = cell(T,N);
    k = cell(T,N);
    L = cell(T,N);
    l = cell(T,N);
    % extract dimensions of system and controls
    n = size(F{1},1);
    all_m = 0;
    for i = 1:N
        m{i} = size(H{1,i},2) - 1 - n;
        start_m{i} =  all_m;
        all_m = all_m + m{i};
    end
    
    
    % Initialize
    for i = 1:N
        G{i} = H{T+1,i}; % constraint-to-go 
        P{i} = Q{T+1,i}; % cost-to-go
    end
    for t = T:-1:1
        Ft = F{t};
        dyn = [1 zeros(1,n+all_m);
               Ft];
        for i = 1:N
            Ht{i} = H{t,i};
            Qt{i} = Q{t,i};
        end

        all_lg = 0;
        all_lh = 0;
        for i = 1:N
            lg{i} = size(G{i},1);
            lh{i} = size(Ht{i},1);
            all_lg = all_lg + lg{i};
            all_lh = all_lh + lh{i};
        end
        
        solve_mat = [];
        for i = 1:N
            M{i} = zeros(m{i} + lh{i} + lg{i} + 2*n,lh{i} + lg{i} + 3*n+1+all_m);
            M{i}(1:m{i},1:1+n+all_m) = Qt{i}(1+n+start_m{i}+1:1+n+start_m{i}+m{i},:);
            M{i}(1:m{i},1+n+all_m+n+1:1+n+all_m+n+n) = Ft(:,1+n+start_m{i}+1:1+n+start_m{i}+m{i})';
            M{i}(1:m{i},1+n+all_m+n+n+1:1+n+all_m+n+n+lh{i}) = Ht{i}(:,1+n+1:end)';
            
            M{i}(m{i}+1:m{i}+n,1) = P{i}(2:end,1);
            M{i}(m{i}+1:m{i}+n,1+n+all_m+1:1+2*n+all_m) = P{i}(2:end,2:end);
            M{i}(m{i}+1:m{i}+n,1+n+all_m+n+1:1+n+all_m+n+n) = -eye(n);
            M{i}(m{i}+1:m{i}+n,1+n+all_m+n+n+lh{i}+1:1+n+all_m+n+n+lh{i}+lg{i}) = G{i}(:,2:end)';
            
            ind = m{i}+n;
            M{i}(ind+1:ind+n,1:1+n+all_m) = Ft;
            M{i}(ind+1:ind+n,1+n+all_m+1:1+2*n+all_m) = -eye(n);
            
            ind = ind+n;
            M{i}(ind+1:ind+lh{i},1:1+n) = Ht{i}(:,1:1+n);
            M{i}(ind+1:ind+lh{i},1+n+start_m{i}+1:1+n+start_m{i}+m{i}) = Ht{i}(:,1+n+1:end);
            
            ind = ind+lh{i};
            M{i}(ind+1:ind+lg{i},1) = G{i}(:,1);
            M{i}(ind+1:ind+lg{i},1+n+all_m+1:1+n+all_m+n) = G{i}(:,2:end);
            
            mu = [M{i}(:,1+n+start_m{i}+1:1+n+start_m{i}+m{i}) M{i}(:,1+n+all_m+1:end)];
            mxv1 = [M{i}(:,1:1+n+start_m{i}) M{i}(:,1+n+start_m{i}+m{i}+1:1+n+all_m)];
            
            pol{i} = -lsqminnorm(mu,...
                                 mxv1);
                             
            lagpol{i} = pol{i}(ind+1:ind+lg{i},:);
            solve_mat = [solve_mat; -pol{i}(1:m{i},1:n+1+start_m{i}) eye(m{i}) -pol{i}(1:m{i},1+n+start_m{i}+1:end)];
        end
        
        allpol = -solve_mat(:,2+n:end)\solve_mat(:,1:1+n);
        
        
        
%         Kmat = lsqminnorm(M,[-Mc -Mx]);
%         allK = Kmat(1:all_m,:);
%        
        ind = 0;
        for i = 1:N
            k{t,i} = allpol(ind+1:ind+m{i},1);
            K{t,i} = allpol(ind+1:ind+m{i},2:end);
            Lag = lagpol{i}*[eye(1+n); allpol(1:start_m{i},:); allpol(start_m{i}+m{i}+1:end,:)];
            L{t,i} = Lag(:,2:end);
            l{t,i} = Lag(:,1);
            ind = ind+m{i};
        end
        
        
%         
        cl = [eye(1+n);
               allpol];
% 
        for i = 1:N
            P{i} = cl'*(Qt{i}+dyn'*P{i}*dyn)*cl;
            G{i} = [Ht{i}(:,1:n+1) zeros(lh{i},start_m{i}) Ht{i}(:,n+2:end) zeros(lh{i},all_m-m{i}-start_m{i});
                                G{i}*dyn]*cl;
            [~,~,V] = svd(G{i});
            rk = rank(G{i},1e-6);
            G{i} = V(:,1:rk)';
        end
    end
end