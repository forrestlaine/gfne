%% Test EC LQ

function [dX,dU,dL,dM,dP] = solve_ec_lq_game_dpol_shared(F,... % time-indexed dict of dynamics
                                                  H,... % time-indexed dict of player-indexed constraints
                                                  HI,...
                                                  Q,... % time-indexed dict of player-indexed costs
                                                  N,... % number of players
                                                  T,... % number of timesteps
                                                  K,...
                                                  m,...
                                                  full_scope)   
    
    WW = cell(T,N);
    ww = cell(T,N);
    % extract dimensions of system and controls
    n = size(F{1},1);
    all_m = 0;
    for i = 1:N
        start_m{i} = all_m;
        all_m = all_m + m{i};
    end
    
    
    % Form base KKT system
    G = H{T+1};
    GI = HI{T+1};
    Ht = H{T};
    HIt = HI{T};
    for i = 1:N
        P{i} = Q{T+1,i};
        Qt{i} = Q{T,i};
    end
    Ft = F{T};
    
    all_lg = size(G,1);
    all_lh = size(Ht,1);
    sHt{T} = all_lh;
    sHt{T+1} = all_lg;
    for i = 1:N
        Gi{i} = G;
        Hti{i} = Ht;
        for j = 1:size(G,1)
            if ~ismember(i,GI{j})
                Gi{i}(j,:) = zeros(1,n+1);
            end
        end
        for j = 1:size(Ht,1)
            if ~ismember(i,HIt{j})
                Hti{i}(j,:) = zeros(1,n+1+all_m);
            end
        end
    end
    
    M{T} = zeros(all_m + all_lh + (N+1)*n + all_lg);
    Mx{T} = zeros(size(M{T},1),n);
    Mc{T} = zeros(size(M{T},1),1);
    
    ind = 0;
    for i = 1:N
        M{T}(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+n+m{i}+1,2+n:end);
        M{T}(ind+1:ind+m{i},all_m+1:all_m+all_lh) = Hti{i}(:,2+n+start_m{i}:1+n+start_m{i}+m{i})';
        M{T}(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind:1+n+ind+m{i})';
        Mx{T}(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+n+m{i}+1,2:1+n);
        Mc{T}(ind+1:ind+m{i}) = Qt{i}(ind+n+2:ind+n+m{i}+1,1);
        ind = ind+m{i};
    end
    
    M{T}(ind+1:ind+all_lh,1:all_m) = Ht(:,2+n:end);
    Mx{T}(ind+1:ind+all_lh,:) = Ht(:,2:1+n);
    Mc{T}(ind+1:ind+all_lh) = Ht(:,1);
    ind = ind+all_lh;
    
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
        M{T}(ind+1:ind+n,ind2+n+1:ind2+n+all_lg) = Gi{i}(:,2:end)';
        Mc{T}(ind+1:ind+n) = P{i}(2:n+1,1);
        ind = ind+n;
    end
    
    M{T}(ind+1:ind+all_lg,ind2+1:ind2+n) = G(:,2:end);
    Mc{T}(ind+1:ind+all_lg) = G(:,1);
    ind = 0;
    
    % Now cycle through remaining steps
    for t = T-1:-1:1
        
        G = H{t+1};
        GI = HI{t+1};
        Ht = H{t};
        HIt = HI{t};
        
        for i = 1:N
            P{i} = Q{t+1,i};
            Qt{i} = Q{t,i};
        end
        Ft = F{t};
        Ftn = F{t+1};

        all_lg = size(G,1);
        all_lh = size(Ht,1);
        sHt{t} = all_lh;
        for i = 1:N
            Gi{i} = G;
            Hti{i} = Ht;
            for j = 1:size(G,1)
                if ~ismember(i,GI{j})
                    Gi{i}(j,:) = zeros(1,n+1+all_m);
                end
            end
            for j = 1:size(Ht,1)
                if ~ismember(i,HIt{j})
                    Hti{i}(j,:) = zeros(1,n+1+all_m);
                end
            end
        end

        top_block = zeros(all_m+all_lh+n+ N*(n+all_m)-all_m,all_m+all_lh+N*n+(N-1)*all_m+n+size(M{t+1},2));
        bleft_block = zeros(size(M{t+1},1), all_m+all_lh+N*n+(N-1)*all_m+n);
        bleft_block(:,end-n+1:end) = Mx{t+1};
        Mx{t} = zeros(size(top_block,1)+size(bleft_block,1),n);
        top_vec = zeros(all_m+all_lh+n+ N*(n+all_m)-all_m,1);
        
        ind = 0;
        for i = 1:N
            top_block(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+n+m{i}+1,2+n:end);
            top_block(ind+1:ind+m{i},all_m+1:all_m+all_lh) = Hti{i}(:,2+n+start_m{i}:1+n+start_m{i}+m{i})';
            top_block(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind:1+n+ind+m{i})';
            Mx{t}(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+n+m{i}+1,2:1+n);
            top_vec(ind+1:ind+m{i}) = Qt{i}(ind+n+2:ind+n+m{i}+1,1);
            ind = ind+m{i};
        end
        ind2 = 0;
        top_block(ind+1:ind+all_lh,1:all_m) = Ht(:,2+n:end);
        Mx{t}(ind+1:ind+all_lh,:) = Ht(:,2:1+n);
        top_vec(ind+1:ind+all_lh) = Ht(:,1);
        ind = ind+all_lh;
        
        top_block(ind+1:ind+n,1:all_m) = Ft(:,2+n:end);
        Mx{t}(ind+1:ind+n,:) = Ft(:,2:1+n);
        top_vec(ind+1:ind+n) = Ft(:,1);
        top_block(ind+1:ind+n,all_m+all_lh+n*N+(N-1)*all_m+1:all_m+all_lh+n*N+n+(N-1)*all_m) = -eye(n);
        ind = ind+n;
    
        ind2 = all_m+all_lh+N*n+(N-1)*all_m+n;
        ind3 = 0;
        ind4 = 0;
        ind5 = 0;
        for i = 1:N
            KK = [];
            for j = 1:N
                if j ~= i
                    KK = [KK; K{t+1,j}];
                end
            end
            top_block(ind+1:ind+(n+all_m-m{i}),ind2-n+1:ind2+all_m) = [P{i}(2:1+n+ind3,2:end);P{i}(1+n+ind3+m{i}+1:end,2:end)];
            top_vec(ind+1:ind+(n+all_m-m{i})) = [P{i}(2:1+n+ind3,1);P{i}(1+n+ind3+m{i}+1:end,1)];
            top_block(ind+1:ind+n,ind2-n-(N-1)*all_m-N*n+1+(i-1)*n:ind2-n-(N-1)*all_m-N*n+i*n) = -eye(n);

            top_block(ind+1:ind+n,ind2-n-(N-1)*all_m+1+ind4:ind2-n-(N-1)*all_m+all_m-m{i}+ind4) = KK';
            top_block(ind+n+1:ind+n+all_m-m{i},ind2-n-(N-1)*all_m+1+ind4:ind2-n-(N-1)*all_m+all_m-m{i}+ind4) = -eye(all_m-m{i});
            
            Gio = Gi{i}(:,2+n:end);
            Gio(:,start_m{i}+1:start_m{i}+m{i}) = [];
            top_block(ind+n+1:ind+n+all_m-m{i},ind2+all_m+1:ind2+all_m+all_lg) = Gio';
            top_block(ind+1:ind+n,ind2+all_m+1:ind2+all_m+all_lg) = Gi{i}(:,2:1+n)';
            if full_scope{i}
                top_block(ind+1:ind+n,ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = Ftn(:,2:n+1)';
            end
            top_block(ind+n+1:ind+(n+all_m-m{i}),ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = [Ftn(:,n+2:n+1+ind3)';Ftn(:,n+2+m{i}+ind3:end)'];
            ind = ind+n+all_m-m{i};
            ind3 = ind3+m{i};
            ind4 = ind4+all_m-m{i};
        end
        M{t} = sparse([top_block; bleft_block M{t+1}]);
        Mc{t} = sparse([top_vec; Mc{t+1}]);
        Mx{t} = sparse(Mx{t});

        ind = 0;
        if t == 1
            if condest(M{t}) > 1/eps
                disp('Rank degenerate in dpol');
            end
            Kmat = M{t}\[-Mc{t} -Mx{t}];
        end
        
%         if t == 1
%             Kmat = lsqminnorm(M{t},);
%         end

%         for i = 1:N
%             ww{t,i} = Kmat(ind+1:ind+m{i},1);
%             WW{t,i} = Kmat(ind+1:ind+m{i},2:end);
%             ind = ind+m{i};
%         end
    end
    
    % Extract solution
    
    sol = Kmat;
%     [dX,dU,dL,dM,dP]

    dX = cell(T+1,1);
    dU = cell(T,N);
    dL = cell(T,N);
    dM = cell(T+1,N);
    dP = cell(T-1,N);
    
    dX{1} = eye(n);
    ind = 0;
    for t = 1:T
        for i = 1:N
           dU{t,i} = sol(ind+1:ind+m{i},:);
           ind = ind+m{i};
        end
        dM{t} = -sol(ind+1:ind+sHt{t},:);
        ind = ind+sHt{t};
        
        for i = 1:N
            dL{t,i} = sol(ind+1:ind+n,:);
            ind = ind+n;
        end
        if t < T
            for i = 1:N
                dP{t,i} = sol(ind+1:ind+all_m-m{i},:);
                ind = ind+all_m-m{i};
            end
        end
        dX{t+1} = sol(ind+1:ind+n,:);
        ind = ind+n;
    end
    dM{T+1} = -sol(ind+1:ind+sHt{T+1},:);
    
end