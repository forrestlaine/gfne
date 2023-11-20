function [x,u,working_set,G_active,flag] = find_feasible_initial_point_lp(x,u,F,H,HI,GG,GI,N,n,m,T,working_set)
    
    all_m = 0;
    
    for i = 1:N
        start_m{i} = all_m;
        all_m = all_m + m{i};
    end
    
    all_p = n*(T) + all_m*T;
    cur_p = [];
    opt = [];
    
    % This might not yield a feasible solution but seems like a good
    % huerustic
    for t = 1:T+1
        

        
        vec = [1;x{t}];
        if t < T+1
            for i = 1:N
                vec = [vec;u{t,i}];
            end
            Gi{t} = zeros(0,1+n+all_m);
            Ge{t} = zeros(0,1+n+all_m);
        else
            Gi{t} = zeros(0,1+n);
            Ge{t} = zeros(0,1+n);
        end
        for j = 1:numel(GG{t})
            if working_set{t}{j}{1} > 0
                k{t,j} = working_set{t}{j}{2};
                Ge{t} = [Ge{t}; GG{t}{j}(k{t,j},:)];
            else
                vals = GG{t}{j}*vec;
                if all(vals < -1e-6)
                    opt = union(opt, GI{t}{j});
                end
                [~,k{t,j}] = max(vals);
                Gi{t} = [Gi{t};GG{t}{j}(k{t,j},:)];
            end
        end
        vals_e = H{t}*vec;
        for j = 1:numel(vals_e)
            if abs(vals_e(j)) > 1e-6
                opt = union(opt, HI{t}{j});
            end
        end
    end
    opt = 1:N;
    dim_not_opt = 0;
    for i = 1:N
        if ~ismember(i,opt)
            dim_not_opt = dim_not_opt + m{i};
        end
    end
    Axeq = [];
    Axineq = [];
    Azeq = [];
    Azineq = [];
    beq = [];
    bineq = [];
    
    ind_x = 0;
    dim_z = 0;
    err = 0;

    for t = 1:T
        
        allu = [];
        for i = 1:N
            allu = [allu; u{t,i}];
        end
        cur_p = [cur_p; allu; x{t+1}];
        nH = size(H{t},1);
        nGe = size(Ge{t},1);
        nGi = size(Gi{t},1);
        block = zeros(n+nH+nGe+dim_not_opt,all_p);
        if t == 1

            block(1:n,ind_x+1:ind_x+all_m) = -F{t}(:,2+n:end);
            block(1:n,ind_x + all_m + 1:ind_x+all_m+n) = eye(n);
            block(n+1:n+nH,ind_x+1:ind_x+all_m) = H{t}(:,2+n:end);
            b = [-F{t}(:,1:1+n)*[1;x{1}]; H{t}(:,1:1+n)*[1;x{1}]];
            d = 0;
            for i = 1:N
                if ~ismember(i,opt)
                    block(n+nH+d+1:n+nH+d+m{i},ind_x+start_m{i}+1:ind_x+start_m{i}+m{i}) = eye(m{i});
                    d = d+m{i};
                    b = [b; -u{t,i}];
                end
            end
            block(n+nH+dim_not_opt+1:end,ind_x+1:ind_x+all_m) = Ge{t}(:,2+n:end);
            b = [b; Ge{t}(:,1:1+n)*[1;x{1}]];
        else
            block(1:n,ind_x+1:ind_x+n+all_m) = -F{t}(:,2:end);
            block(1:n, ind_x+n+all_m+1:ind_x+2*n+all_m) = eye(n);
            block(n+1:n+nH,ind_x+1:ind_x+all_m+n) = H{t}(:,2:end);
            b = [-F{t}(:,1); H{t}(:,1)];
            d = 0;
            for i = 1:N
                if ~ismember(i,opt)
                    block(n+nH+1+d:n+nH+d+m{i},ind_x+n+start_m{i}+1:ind_x+n+start_m{i}+m{i}) = eye(m{i});
                    d = d+m{i};
                    b = [b; -u{t,i}];
                end
            end
            block(n+nH+dim_not_opt+1:end,ind_x+1:ind_x+n+all_m) = Ge{t}(:,2:end);
            b = [b; Ge{t}(:,1)];
        end

        Axeq = sparse([Axeq;block]);
        beq = [beq; b];
        
        block = zeros(nGi,all_p);
        if t == 1
            block(:,ind_x+1:ind_x+all_m) = Gi{t}(:,2+n:end);
            b = Gi{t}(:,1:1+n)*[1;x{1}];
            ind_x = ind_x + all_m;
        else
            block(:,ind_x+1:ind_x+n+all_m) = Gi{t}(:,2:end);
            b = Gi{t}(:,1);
            ind_x = ind_x + n+all_m;
        end
        Axineq = sparse([Axineq; block]);
        bineq = [bineq; b];
        
        
    end
   
    
    nH = size(H{T+1},1);
    nGe = size(Ge{T+1},1);
    nGi = size(Gi{T+1},1);
    
    block = zeros(nH+nGe,all_p);
    block(1:nH,ind_x+1:end) = H{T+1}(:,2:end);
    b = H{T+1}(:,1);
    block(nH+1:end,ind_x+1:end) = Ge{T+1}(:,2:end);
    b = [b; Ge{T+1}(:,1)];
    
    
    Axeq = sparse([Axeq;block]);
    beq = [beq; b];
    
    block = zeros(nGi,all_p);
    block(:,ind_x+1:end) = Gi{T+1}(:,2:end);
    b = Gi{T+1}(:,1);
    Axineq = sparse([Axineq; block]);
    bineq = [bineq; b];    
    
    H = speye(all_p);
    f = -cur_p;
    
    options = optimoptions('quadprog','Display','off');
    [SOL,FVAL,flag] = quadprog(H,f,-Axineq,bineq,-Axeq,beq,[],[],-f,options);
    
    if flag ~= 1
        disp('Did not find a feasible solution');
        working_set = {};
        G_active = {};
        flag = 0;
    else
        flag = 1;
        primals = SOL(1:all_p);
        ind = 0;
        for t = 1:T
            for i = 1:N
                u{t,i} = SOL(ind+1:ind+m{i});
                ind = ind+m{i};
            end
            x{t+1} = SOL(ind+1:ind+n);
            ind = ind+n;
        end
        
        for t = 1:T+1
            G_active{t} = [];
            working_set{t} = {};
            vec = [1;x{t}];
            if t < T+1
                for i = 1:N
                    vec = [vec; u{t,i}];
                end
            end

            for j = 1:numel(GG{t}) 
%                 if false
                vals = GG{t}{j}*vec;
                if all(vals < 1e-6) && t > 1
                    G_active{t} = [G_active{t}; GG{t}{j}(k{t,j},:)];
                    working_set{t}{j}{1} = size(G_active{t},1);
                    working_set{t}{j}{2} = k{t,j};
                    working_set{t}{j}{3} = GI{t}{j};
                else
                    working_set{t}{j}{1} = 0;
                    working_set{t}{j}{2} = 0;
                    working_set{t}{j}{3} = [];
                end
            end
        end
    end
                
            

   
     
end