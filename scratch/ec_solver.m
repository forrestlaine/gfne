%% Iterative EC solver

function [residuals, xval] = ec_solver(f,... % dynamics function
                         h,... % time-indexed dict of player-indexed constraint functions
                         g,... % time-indexed dict of player-indexed inequality constraint functions
                         l,... % time-indexed dict of player-indexed cost functions
                         n,... % dimension of shared state
                         m,... % player-indexed dict of player input dimensions
                         N,... % number of players
                         T,... % time steps
                         z0)   % initial state
    open_loop = false;
    beta = 0.5;
    failed_linesearches = 0;
    import casadi.*
    %% initialize variables and jacobian/hessian functions
    
    shared_constraints = [];

    all_m = 0;
    for i = 1:N
        lagrangian{i} = 0;
        all_m = all_m + m{i};
        control_vars{i} = [];
        other_control_vars{i} = [];
        mu_vars{i} = [];
        lam_vars{i} = [];
        gam_vars{i} = [];
        slack_vars{i} = [];
    end
    
    x{1} = MX.sym(['x_' num2str(1)], n);
    all_vars = [x{1}];
    state_vars = [];
    
    for t = 1:T
        x{t+1} = MX.sym(['x_' num2str(t+1)], n);
        states_controls = [x{t}];
        
        all_vars = [all_vars; x{t+1}];
        state_vars = [state_vars; x{t+1}];

        
        % Constraints
        for i = 1:N
            u{t,i} = MX.sym(['u_' num2str(t) '_' num2str(i)], m{i});
            lam{t,i} = MX.sym(['lam_' num2str(t) '_' num2str(i)], n);
            states_controls = [states_controls; u{t,i}];
            state_single_control = [x{t}; u{t,i}];
            constraint{t,i} = h{t,i}(state_single_control);
            ineq_constraint{t,i} = g{t,i}(state_single_control);
            mu{t,i} = MX.sym(['mu_' num2str(t) '_' num2str(i)], size(constraint{t,i},1));
            gam{t,i} = MX.sym(['gam_' num2str(t) '_' num2str(i)], size(ineq_constraint{t,i},1));
            slacks{t,i} = MX.sym(['slacks_' num2str(t) '_' num2str(i)], size(ineq_constraint{t,i},1));
            
            jac_constraint{t,i} = jacobian(constraint{t,i},state_single_control);
            jac_ineq_constraint{t,i} = jacobian(ineq_constraint{t,i},state_single_control);
                        
            eval_constraint{t,i} = Function('constraint', {state_single_control},{constraint{t,i}});
            eval_jac_constraint{t,i} = Function('jac_constraint', {state_single_control},{jac_constraint{t,i}});
            eval_ineq_constraint{t,i} = Function('ineq_constraint', {state_single_control},{ineq_constraint{t,i}});
            eval_ineq_jac_constraint{t,i} = Function('jac_ineq_constraint', {state_single_control},{jac_ineq_constraint{t,i}});
            
            all_vars = [all_vars; u{t,i}; mu{t,i}; gam{t,i}; slacks{t,i}; lam{t,i}];
            lam_vars{i} = [lam_vars{i}; lam{t,i}];
            control_vars{i} = [control_vars{i}; u{t,i}];
            mu_vars{i} = [mu_vars{i}; mu{t,i}];
            gam_vars{i} = [gam_vars{i}; gam{t,i}];
            slack_vars{i} = [slack_vars{i}; slacks{t,i}];
        end
        
        if t > 1
            for i = 1:N
                for j = 1:N
                    if j ~= i
                        other_control_vars{i} = [other_control_vars{i}; u{t,j}];
                    end
                end
            end
        end
        
        % Dynamics
        pred{t} = f(states_controls);
        jac_pred{t} = jacobian(pred{t},states_controls);
        shared_constraints = [shared_constraints; pred{t}-x{t+1}];
        
        eval_pred{t} = Function('dynamics', {states_controls}, {pred{t}});
        eval_jac_pred{t} = Function('dynamics_jac', {states_controls}, {jac_pred{t}});
        
        % Costs
        for i = 1:N
            cost{t,i} = l{t,i}(states_controls);
            grad_cost{t,i} = gradient(cost{t,i},states_controls);
            hess_cost{t,i} = jacobian(grad_cost{t,i},states_controls);
            
            full_jac_constraint{t,i} = jacobian(constraint{t,i}, states_controls);
            full_jac_ineq_constraint{t,i} = jacobian(ineq_constraint{t,i}, states_controls);
            
            hess_constraint_mult{t,i} = jacobian(full_jac_constraint{t,i}'*mu{t,i},states_controls);
            hess_ineq_constraint_mult{t,i} = jacobian(full_jac_ineq_constraint{t,i}'*gam{t,i},states_controls);
            hess_dyn_mult{t,i} = jacobian(lam{t,i}'*jac_pred{t},states_controls);
            
            eval_cost{t,i} = Function('cost', {states_controls},{cost{t,i}});
            eval_grad_cost{t,i} = Function('cost', {states_controls},{grad_cost{t,i}});
            eval_hess_cost{t,i} = Function('cost', {states_controls},{hess_cost{t,i}});
            
            eval_full_jac_ineq_constraint{t,i} = Function('full_ineq_jac', {states_controls},{full_jac_ineq_constraint{t,i}});
            eval_hess_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,mu{t,i}},{hess_constraint_mult{t,i}});
            eval_hess_ineq_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,gam{t,i}},{hess_ineq_constraint_mult{t,i}});
            eval_hess_dyn_mult{t,i} = Function('hess_dyn', {states_controls,lam{t,i}},{hess_dyn_mult{t,i}});
          
            lagrangian{i} = lagrangian{i} + cost{t,i} +...
                mu{t,i}'*constraint{t,i} +...
                lam{t,i}'*(pred{t}-x{t+1}) -...
                gam{t,i}'*ineq_constraint{t,i};
        end

    end
    final_state = x{T+1};
    for i = 1:N
        constraint{T+1,i} = h{T+1,i}(final_state);
        ineq_constraint{T+1,i} = g{T+1,i}(final_state);
        mu{T+1,i} = MX.sym(['mu_' num2str(T+1) '_' num2str(i)], size(constraint{T+1,i},1));
        gam{T+1,i} = MX.sym(['gam_' num2str(T+1) '_' num2str(i)], size(ineq_constraint{T+1,i},1));
        slacks{T+1,i} = MX.sym(['slacks_' num2str(T+1) '_' num2str(i)], size(ineq_constraint{T+1,i},1));
        
        all_vars = [all_vars; mu{T+1,i}; gam{T+1,i}; slacks{T+1,i}];
        mu_vars{i} = [mu_vars{i}; mu{T+1,i}];
        gam_vars{i} = [gam_vars{i}; gam{T+1,i}];
        slack_vars{i} = [slack_vars{i}; slacks{T+1,i}];
        
        cost{T+1,i} = l{T+1,i}(final_state);
        jac_constraint{T+1,i} = jacobian(constraint{T+1,i},final_state);
        jac_ineq_constraint{T+1,i} = jacobian(ineq_constraint{T+1,i},final_state);
        grad_cost{T+1,i} = gradient(cost{T+1,i},final_state);
        hess_cost{T+1,i} = jacobian(grad_cost{T+1,i},final_state);
        
        hess_constraint_mult{T+1,i} = jacobian(jac_constraint{T+1,i}'*mu{T+1,i},final_state);
        hess_ineq_constraint_mult{T+1,i} = jacobian(jac_ineq_constraint{T+1,i}'*gam{T+1,i},final_state);
        
        eval_constraint{T+1,i} = Function('constraint', {final_state},{constraint{T+1,i}});
        eval_ineq_constraint{T+1,i} = Function('constraint', {final_state},{ineq_constraint{T+1,i}});
        eval_jac_constraint{T+1,i} = Function('jac_constraint', {final_state},{jac_constraint{T+1,i}});
        eval_jac_ineq_constraint{T+1,i} = Function('jac_ineq_constraint', {final_state},{jac_ineq_constraint{T+1,i}});
        eval_cost{T+1,i} = Function('cost', {final_state},{cost{T+1,i}});
        eval_grad_cost{T+1,i} = Function('cost', {final_state},{grad_cost{T+1,i}});
        eval_hess_cost{T+1,i} = Function('cost', {final_state},{hess_cost{T+1,i}});
        eval_hess_constraint_mult{T+1,i} = Function('hess_constraint_mult', {final_state,mu{T+1,i}},{hess_constraint_mult{T+1,i}});
        eval_hess_ineq_constraint_mult{T+1,i} = Function('hess_ineq_constraint_mult', {final_state,gam{T+1,i}},{hess_ineq_constraint_mult{T+1,i}});
        
        lagrangian{i} = lagrangian{i} + cost{T+1,i} + mu{T+1,i}'*constraint{T+1,i} - gam{T+1,i}'*ineq_constraint{T+1,i};
    end
    
    
    eval_dynamics = Function('dyn', {all_vars},{shared_constraints});
    for i = 1:N
        constraint_violation{i} = gradient(lagrangian{i},mu_vars{i}); % h(x)
        ineq_constraint_violation{i} = -slack_vars{i}-gradient(lagrangian{i},gam_vars{i}); % g(x) - s
        
        eval_constraint_violation{i} = Function('con', {all_vars}, {constraint_violation{i}});
        eval_ineq_constraint_violation{i} = Function('con', {all_vars}, {ineq_constraint_violation{i}});
        state_optimality{i} = gradient(lagrangian{i},state_vars);
        eval_state_optimality{i} = Function('con', {all_vars}, {state_optimality{i}});
        control_optimality{i} = gradient(lagrangian{i},control_vars{i});
        eval_control_optimality{i} = Function('con', {all_vars}, {control_optimality{i}});
        other_control_optimality{i} = gradient(lagrangian{i},other_control_vars{i});
        eval_other_control_optimality{i} = Function('con', {all_vars}, {other_control_optimality{i}});
    end
    %% initialize solution from z0 using no control
    xval{1} = z0;
    for t = 1:T
        xu = [xval{t}];
        for i = 1:N
            uval{t,i} = zeros(m{i},1);
            lamval{t,i} = zeros(n,1);
            muval{t,i} = zeros(size(h{t,i}([xval{t};uval{t,i}]),1),1);
            gamval{t,i} = ones(size(g{t,i}([xval{t};uval{t,i}]),1),1);
            slackval{t,i} = ones(size(g{t,i}([xval{t};uval{t,i}]),1),1);
            if t<T
                psival{t,i} = zeros(all_m-m{i},1);
            end
            xu = [xu; uval{t,i}];
        end
        xval{t+1} = f(xu);
    end
    
    for i = 1:N
        muval{T+1,i} = zeros(size(h{T+1,i}(xval{T+1}),1),1);
        gamval{T+1,i} = ones(size(g{T+1,i}(xval{T+1}),1),1);
        slackval{T+1,i} = ones(size(g{T+1,i}(xval{T+1}),1),1);
    end
    %% initialize jacobian/hessian terms
    for t = 1:T
        uu = [];
        for i = 1:N
            H{t,i} = [full(eval_constraint{t,i}([xval{t};uval{t,i}])),...
                      full(eval_jac_constraint{t,i}([xval{t};uval{t,i}]))];
            uu = [uu; uval{t,i}];
        end
        xu = [xval{t};uu];
        F{t} = [full(eval_pred{t}(xu))-xval{t+1}, full(eval_jac_pred{t}(xu))];

        for i = 1:N
            grad = full(eval_grad_cost{t,i}(xu));
            hess = full(eval_hess_cost{t,i}(xu)+...
                   eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                   eval_hess_dyn_mult{t,i}(xu,lamval{t,i}));
            Q{t,i} = [0 grad'; 
                    grad hess];
        end
    end
    for i = 1:N
        H{T+1,i} = [full(eval_constraint{T+1,i}(xval{T+1})), full(eval_jac_constraint{T+1,i}(xval{T+1}))];
        hess = eval_hess_cost{T+1,i}(xval{T+1})+ eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i});
        Q{T+1,i} = [0 full(eval_grad_cost{T+1,i}(xval{T+1}))'; 
                      full(eval_grad_cost{T+1,i}(xval{T+1})) full(hess)];
    end
%     %% initialize policies
    [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);

    %% pack initial z0
    zvec = [z0];
    for t = 1:T
        zvec = [zvec; xval{t+1}];
        for i = 1:N
            zvec = [zvec; uval{t,i}; muval{t,i}; lamval{t,i}];
        end
    end
    for i = 1:N
        zvec = [zvec; muval{T+1,i}];
    end
    cur_necessary_conditions = [eval_dynamics(zvec)];
    for i = 1:N
        % Note that this is valid since our initial estimate of the
        % anticipated-intent multipliers are zero. When these values are
        % non-zero, modifications must be made to the "other_control" and
        % "state" optimality conditions.
        cur_necessary_conditions = [cur_necessary_conditions; 
                                    eval_constraint_violation{i}(zvec); 
                                    eval_state_optimality{i}(zvec);
                                    eval_control_optimality{i}(zvec)];
%                                     eval_other_control_optimality{i}(zvec)];
    end
    current_residual = norm(full(cur_necessary_conditions));
    residuals(1) = current_residual; 
    prev_alpha = 1;
    prev_s = 1;
    %% iteratively solve for solution
    for iter = 1:100 % max iterations
        %% extract update / solve for search direction
%         [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
        [dX,dU,dL,dM,dP] = solve_ec_lq_game_d(F,H,Q,N,T,K,zeros(size(z0)));
%         [dX2,dU2,dL2,dM2,dP2,W,w] = solve_ec_lq_game_both(F,H,Q,N,T,K,k,zeros(size(z0)));
        
            
        %% linesearch (update this so we're not packing and upacking data)
        c1 = 0.1;
        
        alpha = min(1,prev_alpha);   
%         alpha = 1;
        start_s = max(1,prev_s);
%         start_s = 1;
        num_linesearch_iters = 10;

        ls_resids(1) = current_residual;
        linesearch_succeeded = false;
        
        for s = start_s:num_linesearch_iters % max linesearch iterations
            % pack candidate values ('c_') and jacobians, hessians, etc.
            c_xval{1} = z0;
            c_xval2{1} = z0;
            %%
            for t = 1:T
                for i = 1:N
                    c_uval{t,i} = uval{t,i}+alpha*dU{t,i};
                    c_lamval{t,i} = lamval{t,i}+alpha*(dL{t,i}-lamval{t,i});
                    c_muval{t,i} = muval{t,i}+alpha*(dM{t,i}-muval{t,i}); 
                    if t<T
                        c_psival{t,i} = psival{t,i}+alpha*(dP{t,i}-psival{t,i});
                    end
                end
                c_xval{t+1} = xval{t+1}+alpha*dX{t};
                uu = [];
                for i = 1:N
                    c_H{t,i} = [full(eval_constraint{t,i}([c_xval{t};c_uval{t,i}])),...
                                full(eval_jac_constraint{t,i}([c_xval{t};c_uval{t,i}]))];
                    uu = [uu; c_uval{t,i}];
                end
                xu = [c_xval{t};uu];
                c_F{t} = full([full(eval_pred{t}(xu))-c_xval{t+1}, full(eval_jac_pred{t}(xu))]);

                for i = 1:N
                    hess = eval_hess_cost{t,i}(xu)+ ...
                           eval_hess_constraint_mult{t,i}(xu, c_muval{t,i})+ ...
                           eval_hess_dyn_mult{t,i}(xu,c_lamval{t,i});
                    hess(1+1:1+n,1+1:1+n) = hess(1+1:1+n,1+1:1+n);
                    

                    c_Q{t,i} = [0 full(eval_grad_cost{t,i}(xu))'; 
                                full(eval_grad_cost{t,i}(xu)) full(hess)];
                end

            end
            for i = 1:N
                c_muval{T+1,i} = muval{T+1,i}+alpha*(dM{T+1,i}-muval{T+1,i});
                hess = eval_hess_cost{T+1,i}(c_xval{T+1}) + eval_hess_constraint_mult{T+1,i}(c_xval{T+1},c_muval{T+1,i});
                c_H{T+1,i} = [full(eval_constraint{T+1,i}(c_xval{T+1})), full(eval_jac_constraint{T+1,i}(c_xval{T+1}))];
                c_Q{T+1,i} = [0 full(eval_grad_cost{T+1,i}(c_xval{T+1}))'; 
                            full(eval_grad_cost{T+1,i}(c_xval{T+1})) full(hess)];
            end

            % Note these policy jac terms are evaluated after each candidate step
%             c_K = W;
%             if ~open_loop
            [c_K, c_k] = solve_ec_lq_game_r(c_F,c_H,c_Q,N,m,T);
%             end

            new_zvec = [z0];
            for t = 1:T
                new_zvec = [new_zvec; 
                            c_xval{t+1}];
                for i = 1:N
                    new_zvec = [new_zvec; 
                                c_uval{t,i}; 
                                c_muval{t,i}; 
                                c_lamval{t,i}];
                end
            end
            for i = 1:N
                new_zvec = [new_zvec; 
                            c_muval{T+1,i}];
            end


            %% account for feedback constraints in necessary conditions
            new_necessary_conditions = [eval_dynamics(new_zvec)];
            all_conds = [eval_dynamics(new_zvec)];
            state_opt_conds = [];
            for i = 1:N
                % necessary conditions for state and other_control
                % optimalities need to be corrected to account for
                % anticipated-intent constraints/multipliers.
                constraint_opt = eval_constraint_violation{i}(new_zvec);
                state_opt = eval_state_optimality{i}(new_zvec);
                control_opt = eval_control_optimality{i}(new_zvec);
                other_control_opt = eval_other_control_optimality{i}(new_zvec);

  
                ind = 0;
                ind2 = 0;
                for t = 1:T-1
                    c_KK{t+1,i} = [];
                    for j = 1:N
                        if j~=i
                            c_KK{t+1,i} = [c_KK{t+1,i}; c_K{t+1,j}];
                        end
                    end
                    state_opt(ind+1:ind+n) = state_opt(ind+1:ind+n) + c_KK{t+1,i}'*c_psival{t,i};
                    other_control_opt(ind2+1:ind2+all_m-m{i}) = other_control_opt(ind2+1:ind2+all_m-m{i}) - c_psival{t,i};
                    ind = ind+n;
                    ind2 = ind2+all_m-m{i};
                end
                
                if s == num_linesearch_iters
                    disp('check');
                end
                
                new_necessary_conditions = [new_necessary_conditions; 
                                            constraint_opt; 
                                            control_opt;
                                            other_control_opt];
                all_conds = [all_conds; state_opt;
                            constraint_opt; 
                            control_opt;
                            other_control_opt];
                state_opt_conds = [state_opt_conds;
                                    state_opt];
            end

            new_current_residual = norm(full(new_necessary_conditions));
            new_state_resid = norm(full(state_opt_conds))
            all_resid = norm(full(all_conds))
            ls_residual = new_current_residual
            ls_resids(s+1) = new_current_residual;

            if (new_current_residual < current_residual)
                if new_current_residual > current_residual
                    plot_on = true;
                else
                    plot_on = false;
                end
                disp('linesearch accepted!');
                current_residual = new_current_residual;
                prev_alpha = alpha/beta;
                prev_s = s-1;
                linesearch_succeeded = true;
                break;
            else
                disp('linesearch rejected!');
                alpha = alpha * beta;
            end
        end

        %% Accept values after step
        if linesearch_succeeded
            H = c_H;
            Q = c_Q;
            F = c_F;
        end
        
%         if false
%             close all;
%             figure; hold on;
%             spot_a = plot(xval{1}(1), xval{1}(2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%             spot_b = plot(xval{1}(5), xval{1}(6),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%             pause(.1);
%             for t = 1:T
%                 spot_a = plot(xval{t+1}(1), xval{t+1}(2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%                 spot_b = plot(xval{t+1}(5), xval{t+1}(6),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%                 axis([-10 10 -10 10]);
%                 pause(.1);
%             end
%             spot_a = plot(c_xval{1}(1), c_xval{1}(2),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'g');
%             spot_b = plot(c_xval{1}(5), c_xval{1}(6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'g');
%             pause(0.1);
%             for t = 1:T
%                 spot_a = plot(c_xval{t+1}(1), c_xval{t+1}(2),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'g');
%                 spot_b = plot(c_xval{t+1}(5), c_xval{t+1}(6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'g');
%                 axis([-10 10 -10 10]);
%                 pause(0.1);
%             end
%         end
        
%         if alpha < 1e-3
%             hessians_updated = 0;
%             total_hess_norm = 0;
%             for t = 2:T
%                 for i = 1:N
%                     yk = (c_KK{t,i}-KK{t,i})'*psival{t-1,i};
%                     sk = alpha*dX{t};
%                     temp = (yk-HK{t,i}*sk);
%                     if abs(sk'*temp) > 0 && abs(sk'*temp) > norm(sk)*norm(temp)*1e-8
%                         hessians_updated = hessians_updated + 1;
%                         update = temp*temp'/(sk'*temp);
%                         HK{t,i} = HK{t,i}+update;
%                         total_hess_norm = norm(update);
%                         Q{t,i}(1+1:1+n,1+1:1+n) = Q{t,i}(1+1:1+n,1+1:1+n)+update;
%                     end
%                 end
%             end
%         end
        if linesearch_succeeded
            for t = 1:T
                for i = 1:N
                    uval{t,i} = c_uval{t,i};
                    lamval{t,i} = c_lamval{t,i};
                    muval{t,i} = c_muval{t,i}; 
                    if t<T
                        psival{t,i} = c_psival{t,i};
                    end
                end
                xval{t+1} = c_xval{t+1};
                
            end
            for i = 1:N
                muval{T+1,i} = c_muval{T+1,i};
            end
            zvec = new_zvec;
            K = c_K;
            k = c_k;
        else
            disp('ls did not converge!');
            break;
%             failed_linesearches = failed_linesearches + 1;
%             prev_alpha = alpha/beta;
%             prev_s = num_linesearch_iters;
%             
%             [K, ~] = solve_ec_lq_game_r(F,H,Q,N,T);
%             new_necessary_conditions = [eval_dynamics(zvec)];
%             for i = 1:N
%                 % necessary conditions for state and other_control
%                 % optimalities need to be corrected to account for
%                 % anticipated-intent constraints/multipliers.
%                 constraint_opt = eval_constraint_violation{i}(zvec);
%                 state_opt = eval_state_optimality{i}(zvec);
%                 control_opt = eval_control_optimality{i}(zvec);
%                 other_control_opt = eval_other_control_optimality{i}(zvec);
% 
%                 if s == num_linesearch_iters
%                     disp('check');
%                 end
%                 ind = 0;
%                 ind2 = 0;
%                 for t = 1:T-1
%                     KK{t+1,i} = [];
%                     for j = 1:N
%                         if j~=i
%                             KK{t+1,i} = [KK{t+1,i}; K{t+1,j}];
%                         end
%                     end
%                     state_opt(ind+1:ind+n) = state_opt(ind+1:ind+n) + KK{t+1,i}'*psival{t,i};
%                     other_control_opt(ind2+1:ind2+all_m-m{i}) = other_control_opt(ind2+1:ind2+all_m-m{i}) - psival{t,i};
%                     ind = ind+n;
%                     ind2 = ind2+all_m-m{i};
%                 end
%                 
%                new_necessary_conditions = [new_necessary_conditions; 
%                                 state_opt;
%                                 constraint_opt; 
%                                 control_opt;
%                                 other_control_opt];
%                                         
%             end
% 
%             current_residual = norm(full(new_necessary_conditions));
        end
        
        %% update solver stats
        residuals(iter+1) = current_residual;
%         if current_residual < 5
%             close all;
%             figure; hold on;
%             xx(1,:) = z0;
%             spot_a = plot(xx(1,1), xx(1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%             spot_b = plot(xx(1,5), xx(1,6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%             pause(0.1);
%             for t = 1:T
%                 xx(t+1,:) = xval{t+1};
%                 spot_a = plot(xx(t+1,1), xx(t+1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%                 spot_b = plot(xx(t+1,5), xx(t+1,6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%                 axis([-10 10 -10 10]);
%                 pause(0.1);
%             end
%             disp('plotted');
%         end
        %% check for convergence
        current_residual
        if current_residual < 1e-5
            state_optimality=norm(full(state_opt_conds))
            break;
        end
        if failed_linesearches > 10
            break;
        end
    end
    
    close all;
%     writerObj = VideoWriter('ec_fb','MPEG-4');
%     writerObj.FrameRate = 10;
%     open(writerObj);
%     axis tight
%     set(gca,'nextplot','replacechildren');
%     set(gcf,'Renderer','zbuffer');
%     hold on;
    
    half_w = .25;
    half_l = 0.5;
    figure; hold on;
    xx(1,:) = z0;
    spot_a = patch('XData',[xx(1,1)-half_w,xx(1,1)-half_w,xx(1,1)+half_w,xx(1,1)+half_w,xx(1,1)-half_w],...
                   'YData',[xx(1,2)-half_l,xx(1,2)+half_l,xx(1,2)+half_l,xx(1,2)-half_l,xx(1,2)-half_l],...
                   'FaceAlpha',1,'FaceColor', 'b');
    spot_b = patch('XData',[xx(1,5)-half_w,xx(1,5)-half_w,xx(1,5)+half_w,xx(1,5)+half_w,xx(1,5)-half_w],...
                   'YData',[xx(1,6)-half_l,xx(1,6)+half_l,xx(1,6)+half_l,xx(1,6)-half_l,xx(1,6)-half_l],...
                   'FaceAlpha',1,'FaceColor', 'r');
    spot_c = patch('XData',[xx(1,9)-half_w,xx(1,9)-half_w,xx(1,9)+half_w,xx(1,9)+half_w,xx(1,9)-half_w],...
                   'YData',[xx(1,10)-half_l,xx(1,10)+half_l,xx(1,10)+half_l,xx(1,10)-half_l,xx(1,10)-half_l],...
                   'FaceAlpha',1,'FaceColor', 'g');
    rotate(spot_a, [0 0 1], rad2deg(xx(1,4))+90, [xx(1,1) xx(1,2) 0]);
    rotate(spot_b, [0 0 1], rad2deg(xx(1,8))+90, [xx(1,5) xx(1,6) 0]);
    rotate(spot_c, [0 0 1], rad2deg(xx(1,12))+90, [xx(1,9) xx(1,10) 0]);
    
%     spot_a = plot(xx(1,1), xx(1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     spot_b = plot(xx(1,5), xx(1,6),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    axis([-10 10 -10 10]);
    pause(0.1);
    delete(spot_a);
    delete(spot_b);
    delete(spot_c);
    for t = 1:T
        xx(t+1,:) = xval{t+1};
        spot_a = patch('XData',[xx(t+1,1)-half_w,xx(t+1,1)-half_w,xx(t+1,1)+half_w,xx(t+1,1)+half_w,xx(t+1,1)-half_w],...
               'YData',[xx(t+1,2)-half_l,xx(t+1,2)+half_l,xx(t+1,2)+half_l,xx(t+1,2)-half_l,xx(t+1,2)-half_l],...
                   'FaceAlpha',1,'FaceColor', 'b');
        spot_b = patch('XData',[xx(t+1,5)-half_w,xx(t+1,5)-half_w,xx(t+1,5)+half_w,xx(t+1,5)+half_w,xx(t+1,5)-half_w],...
               'YData',[xx(t+1,6)-half_l,xx(t+1,6)+half_l,xx(t+1,6)+half_l,xx(t+1,6)-half_l,xx(t+1,6)-half_l],...
                   'FaceAlpha',1,'FaceColor', 'r');
        spot_c = patch('XData',[xx(t+1,9)-half_w,xx(t+1,9)-half_w,xx(t+1,9)+half_w,xx(t+1,9)+half_w,xx(t+1,9)-half_w],...
                   'YData',[xx(t+1,10)-half_l,xx(t+1,10)+half_l,xx(t+1,10)+half_l,xx(t+1,10)-half_l,xx(t+1,10)-half_l],...
                   'FaceAlpha',1,'FaceColor', 'g');
        rotate(spot_a, [0 0 1], rad2deg(xx(t+1,4))+90, [xx(t+1,1) xx(t+1,2) 0]);
        rotate(spot_b, [0 0 1], rad2deg(xx(t+1,8))+90, [xx(t+1,5) xx(t+1,6) 0]);
        rotate(spot_c, [0 0 1], rad2deg(xx(t+1,12))+90, [xx(t+1,9) xx(t+1,10) 0]);
%         spot_a = plot(xx(t+1,1), xx(t+1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%         spot_b = plot(xx(t+1,5), xx(t+1,6),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        axis([-10 10 -10 10]);
%         frame = getframe;
%         writeVideo(writerObj,frame);
        pause(0.1);
        delete(spot_a);
        delete(spot_b);
        delete(spot_c);
    end
    disp('plotted');
    
%     xx2(1,:) = z0;
%     spot_a = plot(xx2(1,1), xx2(1,2),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'g');
%     spot_b = plot(xx2(1,5), xx2(1,6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'g');
%     pause(0.1);
%     for t = 1:T
%         xx2(t+1,:) = c_xval2{t+1};
%         spot_a = plot(xx2(t+1,1), xx2(t+1,2),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'g');
%         spot_b = plot(xx2(t+1,5), xx2(t+1,6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'g');
%         axis([-10 10 -10 10]);
%         pause(0.1);
%     end
%     close(writerObj);
    sol = current_residual;
end