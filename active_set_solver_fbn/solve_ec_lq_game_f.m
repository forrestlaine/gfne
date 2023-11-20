%% Test EC LQ
% Who knows wtf this version is doing. It tries to extract
% constraints-to-go from KKT of subgame, but it is hard to assign ownership
% of constraints in this formulation. Leads to some weird trajectories that
% appear to be overconstrained and definitely suboptimal.
function [K, k] = solve_ec_lq_game_f(F,... % time-indexed dict of dynamics
                                     H,... % time-indexed dict of player-indexed constraints
                                     Q,... % time-indexed dict of player-indexed costs
                                     N,... % number of players
                                     T)     
    
    K = cell(T,N);
    k = cell(T,N);
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

        M = zeros(all_m + all_lh + (N+1)*n + all_lg);
        Mx = zeros(size(M,1),n);
        Mc = zeros(size(M,1),1);

        ind = 0;
        ind2 = 0;
        ind3 = 0;
        for i = 1:N
            M(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+2:ind+m{i}+1,2+n:end);
            M(ind+1:ind+m{i},all_m+ind2+1:all_m+ind2+lh{i}) = Ht{i}(:,2+n:end)';
            M(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind3:1+n+ind3+m{i})';
            Mx(ind+1:ind+m{i},:) = Qt{i}(ind+2:ind+m{i}+1,2:1+n);
            Mc(ind+1:ind+m{i}) = Qt{i}(ind+2:ind+m{i}+1,1);
            ind = ind+m{i};
            ind2 = ind2+lh{i};
            ind3 = ind3+m{i};
        end
        ind2 = 0;
        for i = 1:N
            M(ind+1:ind+lh{i},1+ind2:1+lh{i}) = Ht{i}(:,2+n:end);
            Mx(ind+1:ind+lh{i},:) = Ht{i}(:,2:1+n);
            Mc(ind+1:ind+lh{i}) = Ht{i}(:,1);
            ind = ind+lh{i};
        end
        M(ind+1:ind+n,1:all_m) = Ft(:,2+n:end);
        Mx(ind+1:ind+n,:) = Ft(:,2:1+n);
        Mc(ind+1:ind+n) = Ft(:,1);
        M(ind+1:ind+n,all_m+all_lh+n*N+1:all_m+all_lh+n*N+n) = -eye(n);
        ind = ind+n;
        ind2 = all_m+all_lh+n*N;
        ind3 = 0;
        for i = 1:N
            M(ind+1:ind+n,ind2+1:ind2+n) = P{i}(2:n+1,2:n+1);
            M(ind+1:ind+n,ind2-(N+1-i)*n+1:ind2-(N-i)*n) = -eye(n);
            M(ind+1:ind+n,ind2+n+ind3+1:ind2+n+ind3+lg{i}) = G{i}(:,2:end)';
            Mc(ind+1:ind+n) = P{i}(2:n+1,1);
            ind = ind+n;
            ind3 = ind3+lg{i};
        end

        for i = 1:N
            M(ind+1:ind+lg{i},ind2+1:ind2+n) = G{i}(:,2:end);
            Mc(ind+1:ind+lg{i}) = G{i}(:,1);
            ind = ind+lg{i};
        end
        ind = 0;
        
%         [U, S, V] = svd(M);
%         rk = rank(M,1e-6);
%         U1 = U(:,1:rk);
%         U2 = U(:,rk+1:end);
%         
%         V1 = V(:,1:rk);
%         V2 = V(rk+1,:);
%         
%         S1 = S(1:rk,1:rk);
%         
%         Kmat = V1/S1*U1'*[-Mc -Mx];
%         resid = U2'*[-Mc -Mx];
        
        Kmat = lsqminnorm(M,[-Mc -Mx]);
        allK = Kmat(1:all_m,:);
       
        for i = 1:N
            k{t,i} = Kmat(ind+1:ind+m{i},1);
            K{t,i} = Kmat(ind+1:ind+m{i},2:end);
            ind = ind+m{i};
        end
        
        pol = [eye(1+n);
               allK];

        for i = 1:N
            P{i} = pol'*(Qt{i}+dyn'*P{i}*dyn)*pol;
            G{i} = [Ht{i}(:,1:n+1) zeros(lh{i},start_m{i}) Ht{i}(:,n+2:end) zeros(lh{i},all_m-m{i}-start_m{i});
                                G{i}*dyn]*pol;
            [~,~,V] = svd(G{i});
            rk = rank(G{i},1e-6);
            G{i} = V(:,1:rk)';
        end
    end
end