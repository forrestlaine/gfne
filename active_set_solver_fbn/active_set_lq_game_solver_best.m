function [x,u,lam,mu,gam,psi,working_set,solve_iters,flag] = active_set_lq_game_solver_best(F,H,HI,G,GI,Q,N,T,m,x,u,working_set,params)
    open_loop = params.open_loop;
    full_scope = params.full_scope;
    solve_iters = 0;
    plot_fn = params.plot_fn;
    
    % working_set{t} = 
    % { { index into G_active (if active, 0 else), [owners of constraint] }, ... }
    
%     for t = T+1:-1:1
%         G_active{t} = [];
%         GI_active{t} = {};
%         for j = 1:numel(G{t})
%             constraint_id = working_set{t}{j}{1};
%             if constraint_id > 0 % active
%                 boundary_id = working_set{t}{j}{2};
%                 owners = working_set{t}{j}{3};
%                 G_active{t} = [G_active{t}; G{t}{j}(boundary_id,:)];
%                 working_set{t}{j}{1} = size(G_active{t},1);
%                 GI_active{t}{size(G_active{t},1)} = owners;
%             end
%         end
%     end
    
    x0 = x{1};
    n = size(x0,1);    
    all_m = 0;
    for i = 1:N
        all_m = all_m + m{i};
    end
    
    [xx,uu,working_set,G_active,flag] = find_feasible_initial_point_lp(x,u,F,H,HI,G,GI,N,n,m,T,working_set);
    if ~flag
        disp('Did not find feasible solution, why?');
        x = 0;
        u=0;
        lam=0;
        mu=0;
        gam=0;
        psi=0;
        working_set=0;
        solve_iters=0;
    else
        x = xx;
        u = uu;

        iters = 0;
        dropped_constraint = false;
        fresh_state = true;
        tried_negs = 0;
        while true
            xsave{iters+1} = x;
            if mod(iters,params.print_frequency) == 0
                disp(['Iteration ' num2str(iters)]);

            end

            failed_attempts = 0;
            while true
                for t = 1:T+1
                    GI_active{t} = {};

                    for j = 1:numel(G{t})
                        if working_set{t}{j}{1} > 0
                            GI_active{t}{working_set{t}{j}{1}} = working_set{t}{j}{3};
                        end
                    end
                    if numel(GI_active{t}) ~= size(G_active{t},1)
                        disp('wtf size wrong');
                    end

                end
                lastwarn('', '');
        
      
               
                [xn,un,lam,mu,gam,psi,flag] = solve_ec_lq_game_super_fast_shared(F,H,HI,G_active,GI_active,Q,N,T,x0,m,full_scope,open_loop);
                [warnMsg, warnId] = lastwarn();
%                 if(~isempty(warnId))
%                     disp('intercept');
%                 end
                if ~flag
                    failed_attempts = failed_attempts + 1;
                    if failed_attempts > 1
                        disp('Very bad');
                    end
                    [x,u,working_set,G_active] = find_feasible_initial_point_lp(x,u,F,H,HI,G,GI,N,n,m,T);
                else
                    break;
                end
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
                        for j = 1:numel(G{t})
                            if working_set{t}{j}{1} > 0 
                                mult_ind = working_set{t}{j}{1};
                                mult_val = gam{t}(mult_ind);
                                c_id = working_set{t}{j}{2};
                                if mult_val < -1e-6
                                    num_neg_mults = num_neg_mults+1;
                                    neg_mult_inds = [neg_mult_inds; [t,j,c_id]];
                                    if mult_val < min_mult
                                        min_mult = mult_val;
                                        min_mult_ind = [t,j,c_id];
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
                    neg_ind = mod(tried_negs, numel(p))+1;
                    ii = p(neg_ind);
                    min_mult_ind = neg_mult_inds(ii,:);
                    % Need to drop a constraint
                    t = min_mult_ind(1);
                    j = min_mult_ind(2);

                    mult_ind = working_set{t}{j}{1};
                    owner = working_set{t}{j}{3};
                    working_set{t}{j}{1} = 0;
                    working_set{t}{j}{2} = 0;
                    working_set{t}{j}{3} = [];
                    if mult_ind == 0
                        disp('check');
                    end
                    for j = 1:numel(G{t})
                        if working_set{t}{j}{1} > mult_ind
                           working_set{t}{j}{1} = working_set{t}{j}{1}-1;
                        end
                    end

                    gam{t}(mult_ind) = [];
                    G_active{t}(mult_ind,:) = [];
                    dropped_constraint = true;
                    dropped_constraint_ind = min_mult_ind;
                    dropped_constraint_owner = owner;
                    tried_negs = tried_negs + 1;
                else
                    % Solution found!
                    disp('solution found');
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

                    for j = 1:numel(G{t})
                        if working_set{t}{j}{1} == 0 % inactive
                            cur_vals = G{t}{j}*[1;vec];
                            if all(cur_vals < -1e-8) 
                                disp('Why is constraint active but not active?');
                            end
                            del_vals = G{t}{j}(:,2:end)*dvec;
                            will_miss = false;
                            alpha_to_enter = 0;
                            alpha_to_leave = inf;
                            [~,index] = min(cur_vals);
                            face_at_entrance = index;
                            for k = 1:size(G{t}{j},1)
                                if cur_vals(k) >= -1e-8
                                    if del_vals(k) < 0 % outside but heading in
                                        alpha_to_enter_face = -(max(cur_vals(k),0)) / del_vals(k);
                                        if alpha_to_enter_face >= alpha_to_enter
                                            alpha_to_enter = alpha_to_enter_face;
                                            face_at_entrance = k;
                                        end
                                    else
                                        will_miss = true;
                                    end
                                else
                                    if del_vals(k) > 0 % inside but leaving
                                        alpha_to_leave = min(alpha_to_leave, -cur_vals(k) / del_vals(k));
                                    end
                                end
                            end
                            if will_miss
                                continue;
                            else
                                if alpha_to_enter < alpha_to_leave
                                    if alpha_to_enter < 1
                                        if face_at_entrance < 1
                                            disp('wtf');
                                        end
                                        violated_inds = [violated_inds; [t,j,face_at_entrance]];
                                    end
                                    if alpha_to_enter < alpha % going to enter constrained region along step
                                        alpha = alpha_to_enter;
                                        min_alpha_ind = [t,j,face_at_entrance];
                                    end
                                end
                            end                           
                        end
                    end
                end
                if alpha < 0
                    disp('VERY BAD! alpha < 0');
                end

                for t = 1:T+1
    %                 dvec = dX{t};
                    cvec = [1;x{t}+alpha*dX{t}];

                    if t <= T
                        for i = 1:N
                            cvec = [cvec;u{t,i}+alpha*dU{t,i}];
                        end
                    end

                    regions_left = [];
                    for j = 1:numel(G{t})
                        if working_set{t}{j}{1} > 0 % active constraint might need dropping
                            if any(G{t}{j}*cvec > 1e-6)
                                % no longer in constraint region
                                cind = working_set{t}{j}{1};
                                regions_left = [regions_left; [t,j,working_set{t}{j}{2}]];
                                if cind > size(gam{t})
                                    disp('why is gam too small?');
                                end
                                gam{t}(cind) = [];
                                G_active{t}(cind,:) = [];
                                working_set{t}{j}{1} = 0;
                                working_set{t}{j}{2} = 0;
                                working_set{t}{j}{3} = [];
                                for k = 1:numel(G{t})
                                    if working_set{t}{k}{1} > cind
                                       working_set{t}{k}{1} = working_set{t}{k}{1}-1;
                                    end
                                end
                                fresh_state = true;
                            end
                        end
                    end

                end

                if alpha < 1
                    t = min_alpha_ind(1);
                    j = min_alpha_ind(2);
                    k = min_alpha_ind(3);

                    if dropped_constraint && all(min_alpha_ind == dropped_constraint_ind)
                        % Just dropped this constraint 
                        if num_neg_mults == 0
                            disp('Num neg mults is 0');
                        end
                        if tried_negs >= num_neg_mults
                            % Have tried all possible drops
                            disp(['Cycle at ', num2str(t), ' ', num2str(j)]);
                            disp(['Tried all possible drops (' num2str(num_neg_mults) ') and none yielded a feasible direction']);
                            neg_mult_inds

%                             disp('Going to just break');
%                             break;

%                             if tried_negs >= 2*num_neg_mults
%                                 disp('Tried everything twice, breaking.');
%                                 break;
%                             end
                            disp(['Number of iters: ', num2str(iters)]);
                       
                            G_active{t} = [G_active{t}; G{t}{j}(k,:)];
                            working_set{t}{j}{1} = size(G_active{t},1);
                            working_set{t}{j}{2} = k;
                            if dropped_constraint_owner == GI{t}{j}
                                working_set{t}{j}{3} = setdiff(1:N,GI{t}{j});
                            else
                                working_set{t}{j}{3} = GI{t}{j};
                            end
                            dropped_constraint = false;
                            fresh_state = true;
                            continue;
                        else
                            G_active{t} = [G_active{t}; G{t}{j}(k,:)];
                            working_set{t}{j}{1} = size(G_active{t},1);
                            working_set{t}{j}{2} = k;
                            working_set{t}{j}{3} = dropped_constraint_owner;
                            dropped_constraint = false;
                        end
                    else
                        % adding different constraint than dropped
                        fresh_state = true;
                    
                        if k > size(G{t}{j},1)
                            disp('what');
                        end
                        G_active{t} = [G_active{t}; G{t}{j}(k,:)];
                        working_set{t}{j}{1} = size(G_active{t},1);
                        working_set{t}{j}{2} = k;
                        working_set{t}{j}{3} = GI{t}{j};
                        dropped_constraint = false;
                    end

                elseif size(regions_left,1) == 0
                    min_mult = inf;
                    min_mult_ind = 0;
                    num_neg_mults = 0;
                    neg_mult_inds = [];
                    for t = 1:T+1
                        for j = 1:numel(G{t})
                            if working_set{t}{j}{1} > 0 
                                mult_ind = working_set{t}{j}{1};
                                c_id = working_set{t}{j}{2};
                                if mult_ind > numel(gam{t})
                                    disp('problem');
                                end
                                mult_val = gam{t}(mult_ind);
                                if mult_val < -1e-6
                                    num_neg_mults = num_neg_mults + 1;
                                    neg_mult_inds = [neg_mult_inds; [t,j,c_id]];
                                    if mult_val < min_mult
                                        min_mult = mult_val;
                                        min_mult_ind = [t,j,c_id];
                                    end
                                end
                            end
                        end
                    end
                    p = randperm(num_neg_mults);
                    tried_negs = 0;
                    fresh_state = false;
                    if min_mult < -1e-6
                        neg_ind = mod(tried_negs, numel(p))+1;
                        ii = p(neg_ind);
                        min_mult_ind = neg_mult_inds(ii,:);
                        % Need to drop a constraint
                        t = min_mult_ind(1);
                        j = min_mult_ind(2);
                        k = working_set{t}{j}{2};
                        mult_ind = working_set{t}{j}{1};
                        owner = working_set{t}{j}{3};
                        working_set{t}{j}{1} = 0;
                        working_set{t}{j}{2} = 0;
                        working_set{t}{j}{3} = [];
                        gam{t}(mult_ind) = [];
                        G_active{t}(mult_ind,:) = [];

                        for j = 1:numel(G{t})
                            if working_set{t}{j}{1} > mult_ind
                                working_set{t}{j}{1} = working_set{t}{j}{1} - 1;
                            end
                        end

                        dropped_constraint = true;
                        dropped_constraint_ind = min_mult_ind;
                        dropped_constraint_owner = owner;
                        tried_negs = tried_negs + 1;
                    else
                        % Make sure assignment inconsistencies are
                        % necessary
                        assignment_mismatch_inds = [];
                        for t=1:T+1
                            for j = 1:numel(G{t})
                                if working_set{t}{j}{1} > 0
                                    if numel(setdiff(working_set{t}{j}{3},GI{t}{j})) ~= 0
                                        assignment_mismatch_inds = [assignment_mismatch_inds; [t,j]];
                                    end
                                end
                            end
                        end
                        if size(assignment_mismatch_inds,1) > 0
                            t = assignment_mismatch_inds(1,1);
                            j = assignment_mismatch_inds(1,2);
                            k = working_set{t}{j}{2};
                            mult_ind = working_set{t}{j}{1};
                            owner = working_set{t}{j}{3};
                            working_set{t}{j}{1} = 0;
                            working_set{t}{j}{2} = 0;
                            working_set{t}{j}{3} = [];
                            gam{t}(mult_ind) = [];
                            G_active{t}(mult_ind,:) = [];
                            dropped_constraint = true;
                            dropped_constraint_ind = [t,j,k];
                            dropped_constraint_owner = owner;
                            continue;
                        else
                        
                        
                            disp('solution found');
                            for t = 1:T+1
                                if t <= T
                                    for i = 1:N
                                        u{t,i} = u{t,i}+alpha*dU{t,i};
                                    end
                                end
                                x{t} = x{t} + alpha*dX{t};
                            end
                            break;
                        end
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


        for t=1:T+1
            for j = 1:numel(G{t})
                if working_set{t}{j}{1} > 0
                    if numel(setdiff(working_set{t}{j}{3},GI{t}{j})) ~= 0
                        disp(['Assignment difference at stage ', num2str(t), ', constraint ', num2str(j), '.']);
                    end
                end
            end
        end

        for t = 1:T+1
            for i = 1:N
                num_ineqs = 0;
                gam_out{t} = zeros(numel(working_set{t}),1);
                for j = 1:numel(working_set{t})
                    if working_set{t}{j}{1} > 0
                        gam_out{t}(j) = gam{t}(working_set{t}{j}{1});
                    else
                        gam_out{t}(j) = 0;
                    end
                end
            end
        end
        gam = gam_out;     
    end
end