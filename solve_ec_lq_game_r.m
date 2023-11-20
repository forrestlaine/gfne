%% Test EC LQ
% This version computes the explicit riccati backup, without considering
% the lagrange multipliers. First substitutes in dynamics, then each player
% minimizes norm of their own constraints (may be suboptimal with respect 
% to their objective if there are multiple ways of doing this).
function [K, k] = solve_ec_lq_game_r(F,... % time-indexed dict of dynamics
                                     H,... % time-indexed dict of player-indexed constraints
                                     Q,... % time-indexed dict of player-indexed costs
                                     N,... % number of players
                                     m,... % dict of player control sizes
                                     T)     
    
    K = cell(T,N);
    k = cell(T,N);
    % extract dimensions of system and controls
    n = size(F{1},1);
    all_m = 0;
   
    for i = 1:N
        m_start{i} = all_m;
        all_m = all_m + m{i};
    end
    
    % Initialize costs-to-go, constraints-to-go
    for i = 1:N
        G{i} = H{T+1,i};
        P{i} = Q{T+1,i};
    end
    
    for t = T:-1:1
        Ft = F{t};
        dyn = [1 zeros(1,n+all_m);
               Ft];
        Mall = [];
        tot = 0;
        for i = 1:N
            Qt{i} = Q{t,i} + dyn'*P{i}*dyn;
            nh = size(H{t,i},1);
            Ht{i} = [H{t,i}(:,1:1+n) zeros(nh, m_start{i}) H{t,i}(:,2+n:end) zeros(nh,all_m-m_start{i}-m{i})];
            Ht{i} = [Ht{i}; G{i}*dyn];
            Hut{i} = Ht{i}(:,1+n+m_start{i}+1:1+n+m_start{i}+m{i});
            [U,S,V] = svd(Hut{i});
            rk{i} = rank(Hut{i},1e-6);
            z{i} = m{i}-rk{i};
            n_start{i} = tot;
            tot = tot+z{i};
            Ur = U(:,1:rk{i});
            Un = U(:,rk{i}+1:end);
            Vr{i} = V(:,1:rk{i});
            Vn{i} = V(:,rk{i}+1:end);
            Sr = S(1:rk{i},1:rk{i});
            Mt{i} = Ur'*Ht{i};
            Mall = [Mall; Mt{i}];
            Rt{i} = Un'*Ht{i};
        end
        
        Vra = blkdiag(Vr{:});
        Vna = blkdiag(Vn{:});
        Va = [Vra Vna];
        
        mod = [eye(1+n) zeros(1+n,all_m);
               zeros(all_m, 1+n) Va];
        % [1;x;u] = mod*[1;x;z1;z2]
        Mall = Mall*mod;
        
        z1pol = -Mall(:,1+n+1:1+n+size(Vra,2))\[Mall(:,1:1+n) Mall(:,2+n+size(Vra,2):end)];
        % z1 = z1pol * [1;x;z2]
        
        z2mod = mod*[eye(1+n) zeros(1+n, size(z1pol,2)-1-n);
                     z1pol;
                     zeros(size(z1pol,2)-1-n,1+n), eye(size(z1pol,2)-1-n)];
        % [1;x;u] = z2mod*[1;x;z2]
       
        Nall = [];
        for i = 1:N
            Qtz{i} = z2mod'*Qt{i}*z2mod;
            Nall = [Nall; Qtz{i}(1+n+n_start{i}+1:1+n+n_start{i}+z{i},:)];
        end
        z2pol = -Nall(:,1+n+1:end)\Nall(:,1:1+n);
        % z2 = z2pol * [1;x]
        
        upol = Va*[z1pol;
                   zeros(size(z2pol,1),1+n) eye(size(z2pol,1))] * [eye(n+1);z2pol];
        ind = 0;
        for i = 1:N
            k{t,i} = upol(ind+1:ind+m{i},1);
            K{t,i} = upol(ind+1:ind+m{i},2:end);
            ind = ind+m{i};
        end
        
        % u = upol*[1;x]
        upd = [eye(1+n);upol];
        rdim = 0;
        for i = 1:N
            P{i} = upd'*Qt{i}*upd;
            G{i} = Rt{i}*upd;
            rdim = rdim + size(G{i},1);
        end
%       
    end
end