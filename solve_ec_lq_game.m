%% Test EC LQ
% This function computes feedback solution when feedback is designed as
% pseudo inverse of full sub game KKT. Probably not correct since other
% players might be minimzing objectives/constraints of other players in
% infeasible case. Also super slow.
function [K, k] = solve_ec_lq_game(F,... % time-indexed dict of dynamics
                                   H,... % time-indexed dict of player-indexed constraints
                                   Q,... % time-indexed dict of player-indexed costs
                                   N,... % number of players
                                   T)     
    
    K = cell(T,N);
    k = cell(T,N);
%     L = cell(T,N);
%     l = cell(T,N);
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
        lg{i} = size(G{i},1);
        lh{i} = size(Ht{i},1);
        all_lg = all_lg + lg{i};
        all_lh = all_lh + lh{i};
    end
    
    M{T} = zeros(all_m + all_lh + (N+1)*n + all_lg);
    Mx{T} = zeros(size(M{T},1),n);
    Mc{T} = zeros(size(M{T},1),1);
    
    ind = 0;
    ind2 = 0;
    ind3 = 0;
    for i = 1:N
        M{T}(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+m{i}+n+1,2+n:end);
        M{T}(ind+1:ind+m{i},all_m+ind2+1:all_m+ind2+lh{i}) = Ht{i}(:,2+n:end)';
        M{T}(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind3:1+n+ind3+m{i})';
        Mx{T}(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+m{i}+n+1,2:1+n);
        Mc{T}(ind+1:ind+m{i}) = Qt{i}(ind+2+n:ind+m{i}+n+1,1);
        ind = ind+m{i};
        ind2 = ind2+lh{i};
        ind3 = ind3+m{i};
    end
    ind2 = 0;
    for i = 1:N
        M{T}(ind+1:ind+lh{i},1+ind2:1+lh{i}) = Ht{i}(:,2+n:end);
        Mx{T}(ind+1:ind+lh{i},:) = Ht{i}(:,2:1+n);
        Mc{T}(ind+1:ind+lh{i}) = Ht{i}(:,1);
        ind = ind+lh{i};
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
        M{T}(ind+1:ind+n,ind2+n+ind3+1:ind2+n+ind3+lg{i}) = G{i}(:,2:end)';
        Mc{T}(ind+1:ind+n) = P{i}(2:n+1,1);
        ind = ind+n;
        ind3 = ind3+lg{i};
    end
    
    for i = 1:N
        M{T}(ind+1:ind+lg{i},ind2+1:ind2+n) = G{i}(:,2:end);
        Mc{T}(ind+1:ind+lg{i}) = G{i}(:,1);
        ind = ind+lg{i};
    end
    ind = 0;
    ind2 = 0;

    Kmat = lsqminnorm(M{T},[-Mc{T} -Mx{T}]);
    
    MT = M{T};
    MxT = Mx{T};
    McT = Mc{T};

    for i = 1:N
        k{T,i} = Kmat(ind+1:ind+m{i},1);
        K{T,i} = Kmat(ind+1:ind+m{i},2:end);
        ind = ind+m{i};
    end
    
    % Now cycle through remaining steps
    for t = T-1:-1:1
        for i = 1:N
            G{i} = H{t+1,i};
            Ht{i} = H{t,i};
            P{i} = Q{t+1,i};
            Qt{i} = Q{t,i};
        end
        Ft = F{t};

        all_lg = 0;
        all_lh = 0;
        for i = 1:N
            lg{i} = size(G{i},1);
            lh{i} = size(Ht{i},1);
            all_lg = all_lg + lg{i};
            all_lh = all_lh + lh{i};
        end

        top_block = zeros(all_m+all_lh+n+ N*(n+all_m)-all_m,all_m+all_lh+N*n+(N-1)*all_m+n+size(M{t+1},2));
        bleft_block = zeros(size(M{t+1},1), all_m+all_lh+N*n+(N-1)*all_m+n);
        bleft_block(:,end-n+1:end) = Mx{t+1};
        Mx{t} = zeros(size(top_block,1)+size(bleft_block,1),n);
        top_vec = zeros(all_m+all_lh+n+ N*(n+all_m)-all_m,1);
        
        ind = 0;
        ind2 = 0;
        ind3 = 0;
        for i = 1:N
            top_block(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+n+m{i}+1,2+n:end);
            top_block(ind+1:ind+m{i},all_m+ind2+1:all_m+ind2+lh{i}) = Ht{i}(:,2+n:end)';
            top_block(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind3:1+n+ind3+m{i})';
            Mx{t}(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+m{i}+n+1,2:1+n);
            top_vec(ind+1:ind+m{i}) = Qt{i}(ind+n+2:ind+m{i}+n+1,1);
            ind = ind+m{i};
            ind2 = ind2+lh{i};
            ind3 = ind3+m{i};
        end
        ind2 = 0;
        for i = 1:N
            top_block(ind+1:ind+lh{i},1+ind2:1+ind2+m{i}) = Ht{i}(:,2+n:end);
            Mx{t}(ind+1:ind+lh{i},:) = Ht{i}(:,2:1+n);
            top_vec(ind+1:ind+lh{i}) = Ht{i}(:,1);
            ind = ind+lh{i};
            ind2 = ind2+m{i};
        end
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
            top_block(ind+n+1:ind+(n+all_m-m{i}),ind2-n-(N-1)*all_m+1+ind4:ind2-n-(N-1)*all_m+all_m-m{i}+ind4) = -eye(all_m-m{i});
            top_block(ind+1:ind+n,ind2+all_m+ind5+1:ind2+all_m+ind5+lg{i}) = G{i}(:,2:1+n)';
            top_block(ind+1:ind+n,ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = Ft(:,2:n+1)';
            top_block(ind+n+1:ind+(n+all_m-m{i}),ind2+all_m+all_lg+(i-1)*n+1:ind2+all_m+all_lg+i*n) = [Ft(:,n+2:n+1+ind3)';Ft(:,n+2+m{i}+ind3:end)'];
            ind = ind+n+all_m-m{i};
            ind3 = ind3+m{i};
            ind4 = ind4+all_m-m{i};
            ind5 = ind5+lg{i};
        end
        M{t} = sparse([top_block; bleft_block M{t+1}]);
        Mc{t} = sparse([top_vec; Mc{t+1}]);
        Mx{t} = sparse(Mx{t});

        ind = 0;
%         Kmat = lsqminnorm(M{t},[-Mc{t} -Mx{t}]);
%         for i = 1:N
%             k{t,i} = Kmat(ind+1:ind+m{i},1);
%             K{t,i} = Kmat(ind+1:ind+m{i},2:end);
%             ind = ind+m{i};
%         end
    end
end