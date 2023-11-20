%% Test EC LQ

function [dX,dU,dL,dM,dP,dD] = solve_ec_lq_game_dev(F,... % time-indexed dict of dynamics
                                               H,... % time-indexed dict of player-indexed constraints
                                               Q,... % time-indexed dict of player-indexed costs
                                               N,... % number of players
                                               T,... % number of timesteps
                                               K,... % precomputed policies
                                               k,...,
                                               x0)   % initial state
    
    % extract dimensions of system and controls
    n = size(F{1},1);
    all_m = 0;
    for i = 1:N
        m{i} = size(H{1,i},2) - 1 - n;
        all_m = all_m + m{i};
    end
    
    
    % Form base KKT system
    for i = 1:N
        G{i} = H{T+1,i};
        Ht{i} = H{T,i};
        P{i} = Q{T+1,i};
        Qt{i} = Q{T,i};
    end
    Ft = F{T};
    
    all_lg = 0;
    all_lh = 0;
    for i = 1:N
        lg{T,i} = size(G{i},1);
        lh{T,i} = size(Ht{i},1);
        all_lg = all_lg + lg{T,i};
        all_lh = all_lh + lh{T,i};
    end
    
    M{T} = zeros(all_m + all_lh + (N+1)*n + all_lg);
    Mx{T} = zeros(size(M{T},1),n);
    Mc{T} = zeros(size(M{T},1),1);
    
    ind = 0;
    ind2 = 0;
    ind3 = 0;
    for i = 1:N
        M{T}(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+n+m{i}+1,2+n:end);
        M{T}(ind+1:ind+m{i},all_m+ind2+1:all_m+ind2+lh{T,i}) = Ht{i}(:,2+n:end)';
        M{T}(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind3:1+n+ind3+m{i})';
        Mx{T}(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+n+m{i}+1,2:1+n);
        Mc{T}(ind+1:ind+m{i}) = Qt{i}(ind+n+2:ind+n+m{i}+1,1);
        ind = ind+m{i};
        ind2 = ind2+lh{T,i};
        ind3 = ind3+m{i};
    end
    ind2 = 0;
    for i = 1:N
        M{T}(ind+1:ind+lh{T,i},1+ind2:1+lh{T,i}) = Ht{i}(:,2+n:end);
        Mx{T}(ind+1:ind+lh{T,i},:) = Ht{i}(:,2:1+n);
        Mc{T}(ind+1:ind+lh{T,i}) = Ht{i}(:,1);
        ind = ind+lh{T,i};
    end
    M{T}(ind+1:ind+n,1:all_m) = Ft(:,2+n:end);
    Mx{T}(ind+1:ind+n,:) = Ft(:,2:1+n);
    Mc{T}(ind+1:ind+n) = Ft(:,1);
    M{T}(ind+1:ind+n,all_m+all_lh+n*N+1:all_m+all_lh+n*N+n) = -eye(n);
    ind = ind+n;
    ind2 = all_m+all_lh+n*N;
    ind3 = 0;
    for i = 1:N
        M{T}(ind+1:ind+n,ind2+1:ind2+n) = P{i}(2:n+1,2:n+1);
        M{T}(ind+1:ind+n,ind2-(N+1-i)*n+1:ind2-(N-i)*n) = -eye(n);
        M{T}(ind+1:ind+n,ind2+n+ind3+1:ind2+n+ind3+lg{T,i}) = G{i}(:,2:end)';
        Mc{T}(ind+1:ind+n) = P{i}(2:n+1,1);
        ind = ind+n;
        ind3 = ind3+lg{T,i};
    end
    
    for i = 1:N
        M{T}(ind+1:ind+lg{T,i},ind2+1:ind2+n) = G{i}(:,2:end);
        Mc{T}(ind+1:ind+lg{T,i}) = G{i}(:,1);
        ind = ind+lg{T,i};
    end
    ind = 0;
    
    % Now cycle through remaining steps
    for t = T-1:-1:1
        for i = 1:N
            G{i} = H{t+1,i};
            Ht{i} = H{t,i};
            P{i} = Q{t+1,i};
            
            Qt{i} = Q{t,i};
        end
        Ft = F{t};
        Ftn = F{t+1};

        all_lg = 0;
        all_lh = 0;
        for i = 1:N
            lg{t,i} = size(G{i},1);
            lh{t,i} = size(Ht{i},1);
            all_lg = all_lg + lg{t,i};
            all_lh = all_lh + lh{t,i};
        end

        top_block = zeros(all_m+all_lh+n+ N*(n+all_m),all_m+all_lh+N*n+N*all_m+n+size(M{t+1},2));
        bleft_block = zeros(size(M{t+1},1), all_m+all_lh+N*n+N*all_m+n);
        bleft_block(:,end-n+1:end) = Mx{t+1};
        Mx{t} = zeros(size(top_block,1)+size(bleft_block,1),n);
        top_vec = zeros(all_m+all_lh+n+ N*(n+all_m),1);
        
        ind = 0;
        ind2 = 0;
        ind3 = 0;
        for i = 1:N
            top_block(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+n+m{i}+1,2+n:end);
            top_block(ind+1:ind+m{i},all_m+ind2+1:all_m+ind2+lh{t,i}) = Ht{i}(:,2+n:end)';
            top_block(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind3:1+n+ind3+m{i})';
            Mx{t}(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+n+m{i}+1,2:1+n);
            top_vec(ind+1:ind+m{i}) = Qt{i}(ind+n+2:ind+n+m{i}+1,1);
            ind = ind+m{i};
            ind2 = ind2+lh{t,i};
            ind3 = ind3+m{i};
        end
        ind2 = 0;
        for i = 1:N
            top_block(ind+1:ind+lh{t,i},ind2+1:ind2+m{i}) = Ht{i}(:,2+n:end);
            Mx{t}(ind+1:ind+lh{t,i},:) = Ht{i}(:,2:1+n);
            top_vec(ind+1:ind+lh{t,i}) = Ht{i}(:,1);
            ind = ind+lh{t,i};
            ind2 = ind2+m{i};
        end
        top_block(ind+1:ind+n,1:all_m) = Ft(:,2+n:end);
        Mx{t}(ind+1:ind+n,:) = Ft(:,2:1+n);
        top_vec(ind+1:ind+n) = Ft(:,1);
        top_block(ind+1:ind+n,all_m+all_lh+n*N+N*all_m+1:all_m+all_lh+n*N+n+N*all_m) = -eye(n);
        ind = ind+n;
    
    
        ind2 = all_m+all_lh+N*n+N*all_m+n;
        ind3 = 0;
        ind4 = 0;
        ind5 = 0;
        all_KK = [];
        all_kk = [];
        for i = 1:N
            KK = [];
            all_KK = [all_KK; K{t+1,i}];
            all_kk = [all_kk; k{t+1,i}];
            for j = 1:N
                if j ~= i
                    KK = [KK; K{t+1,j}];
                end
            end
            top_block(ind+1:ind+(n+all_m-m{i}),ind2-n+1:ind2+all_m) = [P{i}(2:1+n+ind3,2:end);P{i}(1+n+ind3+m{i}+1:end,2:end)];
            top_vec(ind+1:ind+(n+all_m-m{i})) = [P{i}(2:1+n+ind3,1);P{i}(1+n+ind3+m{i}+1:end,1)];
            top_block(ind+1:ind+n,ind2-n-N*all_m-N*n+1+(i-1)*n:ind2-n-N*all_m-N*n+i*n) = -eye(n);

            top_block(ind+1:ind+n,ind2-n-N*all_m+1+ind4:ind2-n-N*all_m+all_m-m{i}+ind4) = KK';
            top_block(ind+n+1:ind+n+all_m-m{i},ind2-n-N*all_m+1+ind4:ind2-n-N*all_m+all_m-m{i}+ind4) = -eye(all_m-m{i});
            top_block(ind+1:ind+n,ind2+all_m+ind5+1:ind2+all_m+ind5+lg{t,i}) = G{i}(:,2:1+n)';
            top_block(ind+1:ind+n,ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = Ftn(:,2:n+1)';
            top_block(ind+n+1:ind+(n+all_m-m{i}),ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = [Ftn(:,n+2:n+1+ind3)';Ftn(:,n+2+m{i}+ind3:end)'];
            ind = ind+n+all_m-m{i};
            ind3 = ind3+m{i};
            ind4 = ind4+all_m-m{i};
            ind5 = ind5+lg{t,i};
        end
        
        top_block(ind+1:ind+all_m, ind2-n-all_m+1:ind2-n) = eye(all_m);
        top_block(ind+1:ind+all_m, ind2-n+1:ind2) = all_KK;
        top_block(ind+1:ind+all_m, ind2+1:ind2+all_m) = -eye(all_m);
        top_vec(ind+1:ind+all_m) = -all_kk;
        
        M{t} = sparse([top_block; bleft_block M{t+1}]);
        Mc{t} = sparse([top_vec; Mc{t+1}]);
        Mx{t} = sparse(Mx{t});

    end
    
    top_block = zeros(N*all_m+n, n+N*(all_m+n)+size(M{1},2));
    bleft_block = zeros(size(M{1},1), n+N*(all_m+n));
    bleft_block(:,end-n+1:end) = Mx{1};
    top_vec = zeros(N*all_m+n,1);
    
    
    
    %%
    ind = 0;
    ind2 = n+N*(all_m+n);
    ind3 = 0;
    ind4 = 0;
    ind5 = 0;
    all_KK = [];
    all_kk = [];
    
    for i = 1:N
        G{i} = H{1,i};
        P{i} = Q{1,i};
    end
    Ftn = F{1};
    
    all_lg = 0;
    for i = 1:N
        lg{1,i} = size(G{i},1);
        all_lg = all_lg + lg{1,i};
    end
    
    top_block(ind+1:ind+n,ind2-n+1:ind2) = -eye(n);
    top_vec(ind+1:ind+n) = -x0;
    ind = ind+n;
    
    for i = 1:N
        KK = [];
        all_KK = [all_KK; K{1,i}];
        all_kk = [all_kk; k{1,i}];
        for j = 1:N
            if j ~= i
                KK = [KK; K{1,j}];
            end
        end
        top_block(ind+1:ind+(n+all_m-m{i}),ind2-n+1:ind2+all_m) = [P{i}(2:1+n+ind3,2:end);P{i}(1+n+ind3+m{i}+1:end,2:end)];
        top_vec(ind+1:ind+(n+all_m-m{i})) = [P{i}(2:1+n+ind3,1);P{i}(1+n+ind3+m{i}+1:end,1)];
        top_block(ind+1:ind+n,ind2-n-N*all_m-N*n+1+(i-1)*n:ind2-n-N*all_m-N*n+i*n) = -eye(n);

        top_block(ind+1:ind+n,ind2-n-N*all_m+1+ind4:ind2-n-N*all_m+all_m-m{i}+ind4) = KK';
        top_block(ind+n+1:ind+n+all_m-m{i},ind2-n-N*all_m+1+ind4:ind2-n-N*all_m+all_m-m{i}+ind4) = -eye(all_m-m{i});
        top_block(ind+1:ind+n,ind2+all_m+ind5+1:ind2+all_m+ind5+lg{1,i}) = G{i}(:,2:1+n)';
        top_block(ind+1:ind+n,ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = Ftn(:,2:n+1)';
        top_block(ind+n+1:ind+(n+all_m-m{i}),ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = [Ftn(:,n+2:n+1+ind3)';Ftn(:,n+2+m{i}+ind3:end)'];
        ind = ind+n+all_m-m{i};
        ind3 = ind3+m{i};
        ind4 = ind4+all_m-m{i};
        ind5 = ind5+lg{1,i};
    end

    top_block(ind+1:ind+all_m, ind2-n-all_m+1:ind2-n) = eye(all_m);
    top_block(ind+1:ind+all_m, ind2-n+1:ind2) = all_KK;
    top_block(ind+1:ind+all_m, ind2+1:ind2+all_m) = -eye(all_m);
    top_vec(ind+1:ind+all_m) = -all_kk;

    M_g = sparse([top_block; bleft_block M{1}]);
    Mc_g = sparse([top_vec; Mc{1}]);
    
    sol = -M_g\Mc_g;
    
    
    %%
    
    % Extract solution
    
%     [dX,dU,dL,dM,dP,dD]

    dX = cell(T,1);
    dU = cell(T,N);
    dD = cell(T,N);
    dP = cell(T,N);
    dL = cell(T+1,N);
    dM = cell(T+1,N);
    
    ind = 0;
    for t = 1:T
        for i = 1:N
            dL{t,i} = sol(ind+1:ind+n);
            ind = ind+n;
        end
        for i = 1:N
            dP{t,i} = sol(ind+1:ind+all_m-m{i});
            ind = ind+all_m-m{i};
        end
        for i = 1:N
            dD{t,i} = sol(ind+1:ind+m{i});
            ind = ind+m{i};
        end
        dX{t} = sol(ind+1:ind+n);
        ind = ind+n;
        for i = 1:N
           dU{t,i} = sol(ind+1:ind+m{i});
           ind = ind+m{i};
        end
        for i = 1:N
            dM{t,i} = sol(ind+1:ind+lh{t,i});
            ind = ind+lh{t,i};
        end
    end
    for i = 1:N
        dL{T+1,i} = sol(ind+1:ind+n);
        ind = ind+n;
    end
    dX{T+1} = sol(ind+1:ind+n);
    ind = ind+n;
    for i = 1:N
        dM{T+1,i} = sol(ind+1:ind+lg{T,i});
        ind = ind+lg{T,i};
    end
    
end