function [x,u,lam,mu,gam,psi,working_set,solve_iters] = active_set_lq_game_solver_mod(F,H,G,GR,Q,N,T,m,x0,working_set,unlimited,params)
    open_loop = params.open_loop;
    contact = params.contact;
    solve_iters = 0;
    H_active = H;
    for t = T+1:-1:1
        for i = 1:N
            for j = 1:size(working_set{t,i},1)
                if working_set{t,i}(j) > 0
                    H_active{t,i} = [H_active{t,i}; G{t,i}(j,:)];
                    working_set{t,i}(j) = size(H_active{t,i},1);
                end
            end
        end
    end
    
    
    
    % WARNING: ASSUMING THE ZERO TRAJ IS FEASIBLE
%     if open_loop
%         [K,k] = solve_ec_lq_game_r(F,H_active,Q,N,m,T);
% %         for t = 1:T
% %             for i = 1:N
% %                 K{t,i} = zeros(size(K{t,i}));
% %             end
% %         end
%         [x,u,lam,mu,psi] = solve_ec_lq_game_d(F,H_active,Q,N,T,K,x0);
%     else
%         [x,u,lam,mu,psi] = solve_ec_lq_game_super_fast_mod(F,H_active,Q,N,T,x0,m);
%     end
    n = size(x0,1);
    for t = 1:T
        x{t} = zeros(n,1);
        for i = 1:N
            u{t,i} = zeros(m{i},1);
        end
    end
    x{T+1} = zeros(n,1);

    all_m = 0;
    for i = 1:N
        all_m = all_m + m{i};
    end
    
    solve_iters = solve_iters+1;
    prevprev_working_set = working_set;
    prev_working_set = working_set;
    while true % until feasible initial point found. 
        feasible = true;
        for t = T+1:-1:1
            vec = [1;x{t}];
            for i = 1:N
                if t < T+1
                    vec = [vec;u{t,i}];
                end
            end
            for i = 1:N
                for j = 1:size(G{t,i},1)
                    if G{t,i}(j,:)*vec < -1e-8 && feasible
                        feasible = false;
                        if t == 1
                            disp('here');
                        end
                        H_active{t,i} = [H_active{t,i}; G{t,i}(j,:)];  
                        working_set{t,i}(j) = size(H_active{t,i},1);
                    end
                end
            end
        end
        if feasible
            break;
        else
            disp('bad! exit!');
            break;
%             if open_loop
%                 [K,k] = solve_ec_lq_game_r(F,H_active,Q,N,m,T);
%                 [x,u,lam,mu,psi] = solve_ec_lq_game_d(F,H_active,Q,N,T,K,x0);
%             else
%                 [x,u,lam,mu,psi,K_old] = solve_ec_lq_game_super_fast_mod(F,H_active,Q,N,T,x0,m,K_old,K_overwrite);
%             end
%             solve_iters = solve_iters+1;
        end
    end
    iters = 0;
    dropped_constraint = false;
    fresh_state = true;
    tried_negs = 0;
    while true
        if iters > 103
            disp('check');
        end
        if open_loop
            [K,k] = solve_ec_lq_game_r(F,H_active,Q,N,m,T);
            [xn,un,lam,mu,psi] = solve_ec_lq_game_d(F,H_active,Q,N,T,K,x0);
        else
            [xn,un,lam,mu,psi,Kk] = solve_ec_lq_game_super_fast_mod(F,H_active,Q,N,T,x0,m);
        end
        solve_iters = solve_iters+1;
        step = 0;
        for t = 1:T+1
            if t <= T
                for i = 1:N
                    dU{t,i} = un{t,i}-u{t,i};
                    step = step + norm(dU{t,i});
                end
            end
            
            dX{t} = xn{t}-x{t};
            step = step + norm(dX{t});
        end
        
        if (step < 1e-6)
            if fresh_state
                min_mult = inf;
                min_mult_ind = 0;
                num_neg_mults = 0;
                neg_mult_inds = [];
                for t = 1:T+1
                    for i = 1:N
                        for j = 1:size(G{t,i},1)
                            if working_set{t,i}(j) > 0 
                                mult_ind = working_set{t,i}(j);
                                mult_val = mu{t,i}(mult_ind);
                                if mult_val < -1e-6
                                    num_neg_mults = num_neg_mults+1;
                                    neg_mult_inds = [neg_mult_inds; [t,i,j]];
                                    if mult_val < min_mult
                                        min_mult = mult_val;
                                        min_mult_ind = [t,i,j];
                                    end
                                end
                            end
                        end
                    end
                end
                p = randperm(num_neg_mults);
                tried_negs = 0;
                fresh_state = false;
            end
            if min_mult < -1e-6
                ii = p(tried_negs+1);
                min_mult_ind = neg_mult_inds(ii,:);
                % Need to drop a constraint
                t = min_mult_ind(1);
                i = min_mult_ind(2);
                j = min_mult_ind(3);
                mult_ind = working_set{t,i}(j);
                working_set{t,i}(j) = 0;
                if mult_ind == 0
                    disp('check');
                end
                mu{t,i}(mult_ind) = [];
                H_active{t,i}(mult_ind,:) = [];
                dropped_constraint = true;
                dropped_constraint_ind = min_mult_ind;
                tried_negs = tried_negs + 1;
            else
                % Solution found!
                break;
            end
        else
            % Need to move, but watch out for constraints
            alpha = 1;
            violated_inds = [];
            for t = 1:T+1
                dvec = dX{t};
                vec = x{t};
                if t <= T
                    for i = 1:N
                        dvec = [dvec;dU{t,i}];
                        vec = [vec; u{t,i}];
                    end
                end
                    
                for i = 1:N
                    for j = 1:size(G{t,i},1)
                        if working_set{t,i}(j) == 0 % inactive
                            if G{t,i}(j,2:end)*dvec < -1e-6 % moving towards boundary
                                % alpha needed to break constraint
                                alpha_to_hit = -(G{t,i}(j,1)+G{t,i}(j,2:end)*vec)/(G{t,i}(j,2:end)*dvec);
                                alpha_to_hit = max(0,alpha_to_hit);
                                if GR{t,i}(j,:)*[1;vec] >= 0 % inside region
                                    if GR{t,i}(j,2:end)*dvec < -1e-6 % leaving region
                                        alpha_to_leave = -(GR{t,i}(j,1)+GR{t,i}(j,2:end)*vec)/(GR{t,i}(j,2:end)*dvec);
                                    else
                                        alpha_to_leave = inf;
                                    end
                                    if alpha_to_hit < alpha_to_leave && alpha_to_hit <= 1
                                        violated_inds = [violated_inds; [t,i,j]];
                                    end
                                    if alpha_to_hit < alpha_to_leave && alpha_to_hit < alpha 
                                        % going to hit constraint in active region
                                        min_alpha_ind = [t,i,j];
                                        alpha = alpha_to_hit;
                                    end
                                else % outside region
                                    if GR{t,i}(j,2:end)*dvec > 1e-6 % entering region
                                        alpha_to_enter = -(GR{t,i}(j,1)+GR{t,i}(j,2:end)*vec)/(GR{t,i}(j,2:end)*dvec);
                                    else
                                        alpha_to_enter = inf;
                                    end
                                    if alpha_to_hit > alpha_to_enter && alpha_to_hit <= 1
                                        violated_inds = [violated_inds; [t,i,j]];
                                    end
                                    if alpha_to_hit > alpha_to_enter && alpha_to_hit < alpha
                                        min_alpha_ind = [t,i,j];
                                        alpha = alpha_to_hit;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            for t = 1:T+1
%                 dvec = dX{t};
                cvec = [1;x{t}+alpha*dX{t}];
                
                if t <= T
                    for i = 1:N
                        cvec = [cvec;u{t,i}+alpha*dU{t,i}];
                    end
                end
                    
                for i = 1:N
                    for j = 1:size(G{t,i},1)
                        if working_set{t,i}(j) > 0 % active constraint might need dropping
                            if GR{t,i}(j,:)*cvec < 0 
                                % no longer in constraint region
                                cind = working_set{t,i}(j);
                                mu{t,i}(cind) = [];
                                H_active{t,i}(cind,:) = [];
                                working_set{t,i}(j) = 0;
                            end
                        end
                    end
                end
            end
            
            if alpha < 1
                t = min_alpha_ind(1);
                i = min_alpha_ind(2);
                j = min_alpha_ind(3);
                    
                if dropped_constraint && all(min_alpha_ind == dropped_constraint_ind)
                    % Just dropped this constraint 
                    if tried_negs == num_neg_mults
                        % Have tried all possible drops
                        disp(['Cycle at ', num2str(t), ' ', num2str(i)]);
                        disp(['Tried all possible drops (' num2str(num_neg_mults) ') and none yielded a feasible direction']);
                        neg_mult_inds
                        
                        if i == 1
                            disp('Not sure how to fix this yet');
                            break;
                        else
                            K_overwrite{t,i} = false;
                            H_active{t,1} = [H_active{t,1}; G{t,i}(j,:)];
                            dropped_constraint = false;
                            fresh_state = true;
                            continue;
                        end
                    end
                else
                    % adding different constraint than dropped
                    fresh_state = true;
                end
                
                H_active{t,i} = [H_active{t,i}; G{t,i}(j,:)];
                working_set{t,i}(j) = size(H_active{t,i},1);
                if size(H_active{t,i}) > m{i}
                    disp('bad!')
                    break;
                end
                dropped_constraint = false;
                if alpha >= 1e-6
                    fresh_state = true;
                end
            else
                min_mult = inf;
                min_mult_ind = 0;
                num_neg_mults = 0;
                neg_mult_inds = [];
                for t = 1:T+1
                    for i = 1:N
                        for j = 1:size(G{t,i},1)
                            if working_set{t,i}(j) > 0 
                                mult_ind = working_set{t,i}(j);
                                if mult_ind > numel(mu{t,i})
                                    disp('problem');
                                end
                                mult_val = mu{t,i}(mult_ind);
                                if mult_val < -1e-6
                                    num_neg_mults = num_neg_mults + 1;
                                    neg_mult_inds = [neg_mult_inds; [t,i,j]];
                                    if mult_val < min_mult
                                        min_mult = mult_val;
                                        min_mult_ind = [t,i,j];
                                    end
                                end
                            end
                        end
                    end
                end
                p = randperm(num_neg_mults);
                tried_negs = 0;
                fresh_state = false;
                if min_mult < -1e-6
                    ii = p(tried_negs+1);
                    min_mult_int = neg_mult_inds(ii,:);
                    % Need to drop a constraint
                    t = min_mult_ind(1);
                    i = min_mult_ind(2);
                    j = min_mult_ind(3);
                    mult_ind = working_set{t,i}(j);
                    working_set{t,i}(j) = 0;
                    mu{t,i}(mult_ind) = [];
                    H_active{t,i}(mult_ind,:) = [];
                    
                    dropped_constraint = true;
                    dropped_constraint_ind = min_mult_ind;
                    tried_negs = tried_negs + 1;
                else
                    % Solution found!
                    break;
                end
            end
            for t = 1:T+1
                if t <= T
                    for i = 1:N
                        u{t,i} = u{t,i}+alpha*dU{t,i};
                    end
                end
                x{t} = x{t} + alpha*dX{t};
            end

        end  
        iters = iters+1;
%         tried_negs
%         if ~unlimited && iters > 800
%             disp('Breaking early!');
%             break;
%         end
    end
    
    
    disp(['Number of iters: ', num2str(iters)]);
    
    for t = 1:T+1
        for i = 1:N
            num_ineqs = 0;
            gam{t,i} = zeros(size(working_set{t,i}));
            for j = 1:size(working_set{t,i},1)
                if working_set{t,i}(j) > 0
                    if working_set{t,i}(j) <= size(mu{t,i},1)
                        num_ineqs = num_ineqs+1;
                        gam{t,i}(j) = mu{t,i}(working_set{t,i}(j));
                    end
                else
                    gam{t,i}(j) = 0;
                end
            end
            mu{t,i} = mu{t,i}(1:end-num_ineqs);
        end
    end
            
          
end