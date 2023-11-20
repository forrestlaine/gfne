%% Iterative EC solver

function [residuals, xval,full_sol] = ec_solver_ip(f,... % dynamics function
                                         h,... % time-indexed dict of player-indexed constraint functions
                                         g,... % time-indexed dict of player-indexed inequality constraint functions
                                         l,... % time-indexed dict of player-indexed cost functions
                                         n,... % dimension of shared state
                                         m,... % player-indexed dict of player input dimensions
                                         N,... % number of players
                                         T,... % time steps
                                         z0,...% initial state
                                         params,...
                                         plot_fn,...
                                         initialization)   
                     
    % Note: the constraints must not depend on decision variables of 
    %       other players, or else this method may not be correct.                 
    import casadi.*
    %% initialize variables and jacobian/hessian functions
    
    shared_constraints = [];

    all_m = 0;
    for i = 1:N
        lagrangian{i} = 0;
        start_m{i} = all_m;
        all_m = all_m + m{i};
        control_vars{i} = [];
        all_control_vars = [];
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
            all_control_vars = [all_control_vars; u{t,i}];
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
        raw_ineq_constraint_violation{i} = -gradient(lagrangian{i},gam_vars{i});
        ineq_constraint_violation{i} = -slack_vars{i}-gradient(lagrangian{i},gam_vars{i}); % g(x) - s
        
        eval_constraint_violation{i} = Function('con', {all_vars}, {constraint_violation{i}});
        eval_ineq_constraint_violation{i} = Function('con', {all_vars}, {ineq_constraint_violation{i}});
        eval_raw_ineq_constraint_violation{i} = Function('con', {all_vars}, {raw_ineq_constraint_violation{i}});
        state_optimality{i} = gradient(lagrangian{i},state_vars);
        eval_state_optimality{i} = Function('con', {all_vars}, {state_optimality{i}});
        control_optimality{i} = gradient(lagrangian{i},control_vars{i});
        full_control_optimality{i} = gradient(lagrangian{i},all_control_vars);
        eval_control_optimality{i} = Function('con', {all_vars}, {control_optimality{i}});
        eval_full_control_optimality{i} = Function('con', {all_vars}, {full_control_optimality{i}});
        other_control_optimality{i} = gradient(lagrangian{i},other_control_vars{i});
        eval_other_control_optimality{i} = Function('con', {all_vars}, {other_control_optimality{i}});
    end
    %% initialize solution from z0 using no control
    for i = 1:N
        tauval{i} = params.tauval_init{i};
    end

    if params.use_initialization % use provided initialization
        xval{1} = z0;
        for t = 1:T-params.init_step
            for i = 1:N
                uval{t,i} = initialization{2}{t,i};
                lamval{t,i} = initialization{3}{t,i};
                muval{t,i} = initialization{4}{t,i};
                slackval{t,i} = max(full(eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),tauval{i});
                gamval{t,i} = tauval{i}./slackval{t,i};
                devval{t,i} = zeros(m{i},1);
                devopt{t,i} = zeros(m{i},1);
                psival{t,i} = initialization{5}{t,i};
                
            end
            xval{t+1} = initialization{1}{t+1};
        end
        for t = T-params.init_step+1:T
            xu = xval{t};
            for i = 1:N
                uval{t,i} = zeros(m{i},1);
                devval{t,i} = zeros(m{i},1);
                devopt{t,i} = zeros(m{i},1);
                lamval{t,i} = zeros(n,1);
                muval{t,i} = zeros(size(h{t,i}([xval{t};uval{t,i}]),1),1);
                slackval{t,i} = max(full(eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),tauval{i});
                gamval{t,i} = tauval{i}./slackval{t,i};
                
                psival{t,i} = zeros(all_m,1);
                xu = [xu; uval{t,i}];
            end
            xval{t+1} = f(xu);
        end
        for i = 1:N
            muval{T+1,i} = zeros(size(h{T+1,i}(xval{T+1}),1),1);
            slackval{T+1,i} = max(full(eval_ineq_constraint{T+1,i}(xval{T+1})),tauval{i});
            gamval{T+1,i} = tauval{i}./slackval{T+1,i};
        end
                
    else % initialize from zero control and zero multipliers
        if params.use_rollout_x0
            xval{1} = params.rollout_x0;
        else
            xval{1} = z0;
        end
        for t = 1:T
            xu = [xval{t}];
            for i = 1:N
                uval{t,i} = zeros(m{i},1);
                devval{t,i} = zeros(m{i},1);
                devopt{t,i} = zeros(m{i},1);
                lamval{t,i} = zeros(n,1);
                muval{t,i} = zeros(size(h{t,i}([xval{t};uval{t,i}]),1),1);
                slackval{t,i} = max(full(eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),tauval{i});
                gamval{t,i} = tauval{i}./slackval{t,i};
                psival{t,i} = zeros(all_m,1);
                xu = [xu; uval{t,i}];
            end
            xval{t+1} = f(xu);
        end
        xval{1} = z0;
    end
    
    for i = 1:N
        muval{T+1,i} = zeros(size(h{T+1,i}(xval{T+1}),1),1);
        slackval{T+1,i} = max(full(eval_ineq_constraint{T+1,i}(xval{T+1})),tauval{i});
        gamval{T+1,i} = tauval{i}./slackval{T+1,i};
    end

    %% Now initialize jacobian/hessian terms for real
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
            full_ineq_jac{t,i} = eval_full_jac_ineq_constraint{t,i}(xu);
            Sig = sparse(diag(gamval{t,i}./slackval{t,i}));
            ineq_const_val{t,i} = eval_ineq_constraint{t,i}([xval{t};uval{t,i}]);
            grad = full(eval_grad_cost{t,i}(xu) +...
                        full_ineq_jac{t,i}'*(Sig*ineq_const_val{t,i} -...
                                        gamval{t,i} -...
                                        tauval{i} * slackval{t,i}.^(-1)));
            hess = eval_hess_cost{t,i}(xu)+...
                   eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                   eval_hess_dyn_mult{t,i}(xu,lamval{t,i})-...
                   eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t,i});
            aug = full_ineq_jac{t,i}'*Sig*full_ineq_jac{t,i};
            hess = full(hess+aug);
            
            Q{t,i} = [0 grad'; 
                      grad hess];
        end
    end
    for i = 1:N
        H{T+1,i} = [full(eval_constraint{T+1,i}(xval{T+1})), full(eval_jac_constraint{T+1,i}(xval{T+1}))];
        
        full_ineq_jac{T+1,i} = eval_jac_ineq_constraint{T+1,i}(xval{T+1});
        Sig = sparse(diag(gamval{T+1,i}./slackval{T+1,i}));
        ineq_const_val{T+1,i} = eval_ineq_constraint{T+1,i}(xval{T+1});
        grad = full(eval_grad_cost{T+1,i}(xval{T+1}) +...
                        full_ineq_jac{T+1,i}'*(Sig*ineq_const_val{T+1,i} -...
                                        gamval{T+1,i} -...
                                        tauval{i} * slackval{T+1,i}.^(-1)));
            
            
        hess = eval_hess_cost{T+1,i}(xval{T+1})+...
            eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i})-...
            eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1,i});
        aug = full_ineq_jac{T+1,i}'*Sig*full_ineq_jac{T+1,i};
        hess = full(hess+aug);
        Q{T+1,i} = [0 grad'; 
                    grad hess];
    end  
    %% Now initialize policies
    [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
    fresh_policies = true;
    fresh_policy_count = 1;
    if params.open_loop || params.wait_for_ol_convergence
        for t = 1:T
            for i = 1:N
                K{t,i} = zeros(size(K{t,i}));
                k{t,i} = zeros(size(k{t,i}));
            end
        end
        fresh_policies = false;
        fresh_policy_count = 0;
    end

    %% pack initial z0
    zvec = [z0];
    for t = 1:T
        zvec = [zvec; xval{t+1}];
        for i = 1:N
            zvec = [zvec; uval{t,i}; muval{t,i}; gamval{t,i}; slackval{t,i}; lamval{t,i}];
        end
    end
    for i = 1:N
        zvec = [zvec; muval{T+1,i}; gamval{T+1,i}; slackval{T+1,i}];
    end

    current_residual = inf;
    current_residual_partial = inf;
    residuals(1) = current_residual; 
    prev_alpha = 1;
    prev_s = 1;
    frac_to_bound = 0.995;
    tauval_decreased = false;
    ol_converged = ~params.wait_for_ol_convergence;
    cl_failed = false;
    restoration = false;
    feedback_refinement = false;
    reg = 0;
    %% iteratively solve for solution
    for iter = 1:100 % max iterations
        %% extract update / solve for search direction
        tic;
%         [dXq,dUq,dLq,dMq,dPq] = solve_ec_lq_game_d(F,H,Q,N,T,K,zeros(size(z0)));
        
        
        for t = 1:T
            for i = 1:N
                kd{t,i} = K{t,i}*xval{t}+k{t,i}+devval{t,i}-uval{t,i};
                devdelt{t,i} = devval{t,i}-devopt{t,i};
            end
        end
%         if feedback_refinement
%             reg = params.dev_penalty;
%         else
%             reg = 0;
%         end
        [dX,dU,dL,dM,dP,dD] = solve_ec_lq_game_regdev(F,H,Q,N,T,K,kd,zeros(size(z0)),reg,devdelt);
        

        
        
        solve_time = toc;
%         [dX2,dU2,dL2,dM2,dP2,dD2] = solve_ec_lq_game_dev(F,H,Q,N,T,K,kd,zeros(size(z0)));
%         [dX2,dU2,dL2,dM2,dP2,W,w] = solve_ec_lq_game_both(F,H,Q,N,T,K,k,zeros(size(z0)));
        for i = 1:N
            alpha_dual{i} = 1;
            alpha_primal{i} = 1;
        end
        tic;
        for t = 1:T+1
            dxu = dX{t};
            if t < T+1
                for j = 1:N
                    dxu = [dxu; dU{t,i}];
                end
            end
                
            for i = 1:N
                Sig = sparse(diag(gamval{t,i}./slackval{t,i}));
                dG{t,i} = full(-Sig*(full_ineq_jac{t,i}*dxu + ineq_const_val{t,i} - slackval{t,i} - tauval{i}*gamval{t,i}.^(-1)));
                dS{t,i} = full(Sig\(tauval{i}*slackval{t,i}.^(-1) - dG{t,i}));
                
                for j = 1:size(dG{t,i})
                    if frac_to_bound*gamval{t,i}(j) + alpha_dual{i}*(dG{t,i}(j)-gamval{t,i}(j)) < 1e-12
                        alpha_dual{i} = (frac_to_bound*gamval{t,i}(j)-1e-12)/(gamval{t,i}(j)-dG{t,i}(j));
                        
                    end
                    if frac_to_bound*slackval{t,i}(j) + alpha_primal{i}*dS{t,i}(j) < 1e-12
                        alpha_primal{i} = (1e-12 -frac_to_bound*slackval{t,i}(j))/dS{t,i}(j);
                    end
                end 
            end
        end
        for i = 1:N
            if alpha_primal{i} < 0 || alpha_dual{i} < 0 || alpha_dual{i} > 1 || alpha_primal{i} > 1
                disp('VERY BAD');
            end
        end
        extract_slack_time = toc;

        %% linesearch (update this so we're not packing and upacking data)
       
        if (iter == 1 && params.crash_start)
            restoration = true;
            for i = 1:N
                alpha_primal{i} = 1;
                alpha_dual{i} = 1;
            end
        end
        if iter <= params.warm_up_iters
            alpha = params.warm_up_alpha;
        else
            alpha = params.alpha_init;
        end
        start_s = 1;
        
        if ~params.independent_updates
            min_alpha_primal = inf;
            min_alpha_dual = inf;
            for i = 1:N
                min_alpha_primal = min(alpha_primal{i},min_alpha_primal);
                min_alpha_dual = min(alpha_dual{i},min_alpha_dual);
            end
            for i = 1:N
                alpha_primal{i} = min_alpha_primal;
                alpha_dual{i} = min_alpha_dual;
            end
            if params.single_alpha
                alpha_both = min(alpha_primal{1},alpha_dual{1});
                for i = 1:N
                    alpha_primal{i} = alpha_both;
                    alpha_dual{i} = alpha_both;
                end
            end
        end
        
        if params.debug_print >= 2
            alpha_primal
            alpha_dual
        end

        ls_resids(1) = current_residual;
        linesearch_succeeded = false;
        if params.open_loop || ~ol_converged
            num_linesearch_iters = params.max_linesearch_iters;
        else
            num_linesearch_iters = params.max_cl_linesearch_iters;
        end
        if cl_failed
            num_linesearch_iters = 1;
            alpha = 0.1;
        end
        tic;
        for s = start_s:num_linesearch_iters % max linesearch iterations
            % pack candidate values ('c_') and jacobians, hessians, etc.
            c_xval{1} = z0;
            %%
            for t = 1:T
                for i = 1:N
                    c_uval{t,i} = uval{t,i}+alpha*alpha_primal{i}*dU{t,i};
                    c_devval{t,i} = devval{t,i}+alpha*alpha_primal{i}*dD{t,i};
                    c_lamval{t,i} = lamval{t,i}+alpha*alpha_dual{i}*(dL{t+1,i}-lamval{t,i});
                    c_muval{t,i} = muval{t,i}+alpha*alpha_dual{i}*(dM{t,i}-muval{t,i}); 
%                     if t<T
                    c_psival{t,i} = psival{t,i}+alpha*alpha_dual{i}*(dP{t,i}-psival{t,i});
%                     end
                    if params.feasible || restoration
                        c_slackval{t,i} = max(full(eval_ineq_constraint{t,i}([c_xval{t};c_uval{t,i}])),params.slack_min);
                    else
                        c_slackval{t,i} = slackval{t,i}+alpha*alpha_primal{i}*dS{t,i};
                    end
                    if restoration 
                        c_gamval{t,i} = params.tauval_tolerance./c_slackval{t,i};
                    else
                        c_gamval{t,i} = gamval{t,i}+alpha*alpha_dual{i}*(dG{t,i}-gamval{t,i});
                    end
                    
                    if any(c_gamval{t,i} < 0) || any(c_slackval{t,i} < 0)
                        c_gamval{t,i} = max(0.01,c_gamval{t,i});
                        c_slackval{t,i} = max(0.01,c_slackval{t,i});
                        disp('Warning! Correcting negative slacks or mults');
                    end
                end
                c_xval{t+1} = xval{t+1}+alpha*alpha_primal{i}*dX{t+1};
                
%                 uu = [];
%                 for i = 1:N
%                     c_H{t,i} = [full(eval_constraint{t,i}([c_xval{t};c_uval{t,i}])),...
%                                 full(eval_jac_constraint{t,i}([c_xval{t};c_uval{t,i}]))];
%                     uu = [uu; c_uval{t,i}];
%                 end
%                 xu = [c_xval{t};uu];
%                 c_F{t} = full([full(eval_pred{t}(xu))-c_xval{t+1}, full(eval_jac_pred{t}(xu))]);

%                 for i = 1:N
%                     c_full_ineq_jac{t,i} = eval_full_jac_ineq_constraint{t,i}(xu);
%                     Sig = sparse(diag(c_gamval{t,i}./c_slackval{t,i}));
%                     c_ineq_const_val{t,i} = eval_ineq_constraint{t,i}([c_xval{t};c_uval{t,i}]);
%                     grad = full(eval_grad_cost{t,i}(xu) +...
%                                 c_full_ineq_jac{t,i}'*(Sig*c_ineq_const_val{t,i} -...
%                                                 c_gamval{t,i} -...
%                                                 tauval{i} * c_slackval{t,i}.^(-1)));
%                     hess = eval_hess_cost{t,i}(xu)+...
%                            eval_hess_constraint_mult{t,i}(xu,c_muval{t,i})+...
%                            eval_hess_dyn_mult{t,i}(xu,lamval{t,i})-...
%                            eval_hess_ineq_constraint_mult{t,i}(xu,c_gamval{t,i});
%                     aug = c_full_ineq_jac{t,i}'*Sig*c_full_ineq_jac{t,i};
%                     hess = full(hess+aug);
% 
%                     c_Q{t,i} = [0 grad'; 
%                                 grad hess];
%                     
%                 end

            end
            for i = 1:N
                c_muval{T+1,i} = muval{T+1,i}+alpha*alpha_dual{i}*(dM{T+1,i}-muval{T+1,i});
                
                if params.feasible || restoration
                    c_slackval{T+1,i} = max(full(eval_ineq_constraint{T+1,i}(c_xval{T+1})),params.slack_min);
                else
                    c_slackval{T+1,i} = slackval{T+1,i}+alpha*alpha_primal{i}*dS{T+1,i};
                end
                if restoration
                    c_gamval{T+1,i} = params.tauval_tolerance./c_slackval{T+1,i};
                else
                    c_gamval{T+1,i} = gamval{T+1,i}+alpha*alpha_dual{i}*(dG{T+1,i}-gamval{T+1,i});
                end
                if any(c_gamval{T+1,i} < 0) || any(c_slackval{T+1,i} < 0)
                    c_gamval{T+1,i} = max(0.01,c_gamval{T+1,i});
                    c_slackval{T+1,i} = max(0.01,c_slackval{T+1,i});
                    disp('Warning! Correcting negative slacks or mults');
                end
                
%                 c_H{T+1,i} = [full(eval_constraint{T+1,i}(c_xval{T+1})), full(eval_jac_constraint{T+1,i}(c_xval{T+1}))];
%                 
%                 c_full_ineq_jac{T+1,i} = eval_jac_ineq_constraint{T+1,i}(c_xval{T+1});
%                 Sig = sparse(diag(c_gamval{T+1,i}./c_slackval{T+1,i}));
%                 c_ineq_const_val{T+1,i} = eval_ineq_constraint{T+1,i}(c_xval{T+1});
%                 grad = full(eval_grad_cost{T+1,i}(c_xval{T+1}) +...
%                                 c_full_ineq_jac{T+1,i}'*(Sig*c_ineq_const_val{T+1,i} -...
%                                                 c_gamval{T+1,i} -...
%                                                 tauval{i} * c_slackval{T+1,i}.^(-1)));
% 
% 
%                 hess = eval_hess_cost{T+1,i}(c_xval{T+1})+...
%                     eval_hess_constraint_mult{T+1,i}(c_xval{T+1},c_muval{T+1,i})-...
%                     eval_hess_ineq_constraint_mult{T+1,i}(c_xval{T+1},c_gamval{T+1,i});
%                 aug = c_full_ineq_jac{T+1,i}'*Sig*c_full_ineq_jac{T+1,i};
%                 hess = full(hess+aug);
%                 c_Q{T+1,i} = [0 grad'; 
%                             grad hess];
                        
%                 hess = eval_hess_cost{T+1,i}(c_xval{T+1}) + eval_hess_constraint_mult{T+1,i}(c_xval{T+1},c_muval{T+1,i});
%                 c_Q{T+1,i} = [0 full(eval_grad_cost{T+1,i}(c_xval{T+1}))'; 
%                             full(eval_grad_cost{T+1,i}(c_xval{T+1})) full(hess)];
            end

            % Note these policy jac terms are evaluated after each candidate step
%             c_K = W;
            if params.open_loop || ~ol_converged 
                c_K = K;
                c_k = k;
            else
                % TODO obvi won't work with removed above
                c_K = K;
                c_k = k;
%                 [c_K, c_k] = solve_ec_lq_game_r(c_F,c_H,c_Q,N,m,T);
            end

            new_zvec = [z0];
            for t = 1:T
                new_zvec = [new_zvec; 
                            c_xval{t+1}];
                for i = 1:N
                    new_zvec = [new_zvec; 
                                c_uval{t,i}; 
                                c_muval{t,i}; 
                                c_gamval{t,i};
                                c_slackval{t,i};
                                c_lamval{t,i}];
                end
            end
            for i = 1:N
                new_zvec = [new_zvec; 
                            c_muval{T+1,i};
                            c_gamval{T+1,i};
                            c_slackval{T+1,i}];
            end


            %% account for feedback constraints in necessary conditions
            dynamic_conds = eval_dynamics(new_zvec);
            all_constraints = [dynamic_conds];
            new_necessary_conditions = [dynamic_conds];
            all_conds = [dynamic_conds];
            state_opt_conds = [];
            for i = 1:N
                % necessary conditions for state and other_control
                % optimalities need to be corrected to account for
                % anticipated-intent constraints/multipliers.
                constraint_opt = eval_constraint_violation{i}(new_zvec);
                ineq_constraint_opt = eval_ineq_constraint_violation{i}(new_zvec);
                raw_ineq_constraint_opt = eval_ineq_constraint_violation{i}(new_zvec);
                state_opt = eval_state_optimality{i}(new_zvec);
                control_opt = eval_control_optimality{i}(new_zvec);
                all_control_opt = eval_full_control_optimality{i}(new_zvec);
                other_control_opt = eval_other_control_optimality{i}(new_zvec);
                complementarity_opt = [];
                policy_opt = [];
                deviation_opt = [];

                ind = 0;
                ind2 = 0;
    
                for t = 1:T
                    c_KK = [];
                    c_kk = [];
                    for j = 1:N
                        c_KK = [c_KK; c_K{t,j}];
                        c_kk = [c_kk; c_k{t,j}];
                    end
                    deviation_opt = [deviation_opt; 
                                     reg*(c_devval{t,i}-devopt{t,i}) + c_psival{t,i}(start_m{i}+1:start_m{i}+m{i})];
%                    
                    if t > 1
                        state_opt(ind+1:ind+n) = state_opt(ind+1:ind+n) + c_KK'*c_psival{t,i};
                        ind = ind+n;
                    end
                    
                    policy_opt = [policy_opt; K{t,i}*c_xval{t}+k{t,i}+c_devval{t,i}-c_uval{t,i}];
                    
                    all_control_opt(ind2+1:ind2+all_m) = all_control_opt(ind2+1:ind2+all_m) - c_psival{t,i};
                    ind2 = ind2+all_m;
                end
                for t= 1:T+1
                    complementarity_opt = [complementarity_opt; c_gamval{t,i}.*c_slackval{t,i}-tauval{i}];
                end
                private_conds{i} = [constraint_opt;
                                    state_opt;
                                    ineq_constraint_opt;
                                    all_control_opt;
                                    complementarity_opt;
                                    deviation_opt;
                                    policy_opt];
                c_private_resids{i,1} = norm(full(constraint_opt));
                c_private_resids{i,2} = norm(full(state_opt));
                c_private_resids{i,3} = norm(full(ineq_constraint_opt));
                c_private_resids{i,4} = norm(min(full(raw_ineq_constraint_opt),0));
                c_private_resids{i,5} = norm(full(all_control_opt));
                c_private_resids{i,6} = norm(full(complementarity_opt));
                c_private_resids{i,7} = norm(full(deviation_opt));
                c_private_resids{i,8} = norm(full(policy_opt));
                
                
                all_constraints = [all_constraints;
                                   min(full(eval_raw_ineq_constraint_violation{i}(new_zvec)),0);
                                   constraint_opt];
                
                new_necessary_conditions = [new_necessary_conditions; 
                                            private_conds{i}];
                all_conds = [all_conds; 
                            state_opt;
                            private_conds{i}];
                        
                private_conds{i} = [private_conds{i}; 
                                    dynamic_conds];
                state_opt_conds = [state_opt_conds;
                                    state_opt];
                new_private_residual{i} = norm(full(private_conds{i}));
                if s == num_linesearch_iters
                    disp('check');
                end
            end
            
            new_current_residual = norm(full(new_necessary_conditions));
            only_constraint_residual = norm(full(all_constraints));
            
%             new_state_resid = norm(full(state_opt_conds))
%             all_resid = norm(full(all_conds))
%             current_residual
%             ls_residual = new_current_residual_ls
            ls_resids(s+1) = new_current_residual;
            if (new_current_residual < current_residual || tauval_decreased )
%             if (new_current_residual < current_residual || tauval_decreased || (fresh_policies && params.converge_on_every_policy) )
                if params.debug_print >= 1
                    disp('linesearch accepted!');
                end
                restoration = false;
                break;
            else
                if params.debug_print >= 1
                    disp('linesearch rejected!');
                end
                alpha = alpha * params.beta;
                if s == num_linesearch_iters
                    failed_linesearch = true;
                    if params.restore_after_failed_ls
                        restoration = true;
                        for i = 1:N
                            tauval{i} = params.tauval_init{i};
                        end
                    end
                end
            end
        end
        ls_time = toc;
%         alpha_primal
%         alpha_dual
        
        current_residual = new_current_residual;
        private_residual = new_private_residual;
        private_fine_resids = c_private_resids;
        if params.debug_print >= 1
            current_residual
        end

%         current_residual_partial = new_current_residual_partial;
%         H = c_H;
%         Q = c_Q;
%         F = c_F;
%         full_ineq_jac = c_full_ineq_jac;
%         ineq_const_val = c_ineq_const_val;
        for i = 1:N
            min_comp{i} = 0;
            comp{i} = 0;
            dim{i} = 0;
        end
        for t = 1:T
            for i = 1:N
                uval{t,i} = c_uval{t,i};
                devval{t,i} = c_devval{t,i};
                lamval{t,i} = c_lamval{t,i};
                muval{t,i} = c_muval{t,i}; 
                if t<T
                    psival{t,i} = c_psival{t,i};
                end
                gamval{t,i} = c_gamval{t,i};
                slackval{t,i} = c_slackval{t,i};
                if size(gamval{t,i},1) > 0
                    comp{i} = comp{i} + gamval{t,i}'*slackval{t,i};
                    min_comp{i} = min(min_comp{i}, min(gamval{t,i}.*slackval{t,i}));
                end
                dim{i} = dim{i} + size(gamval{t,i},1);
                
            end
            xval{t+1} = c_xval{t+1};

        end
        for i = 1:N
            muval{T+1,i} = c_muval{T+1,i};
            gamval{T+1,i} = c_gamval{T+1,i};
            slackval{T+1,i} = c_slackval{T+1,i};
            if size(gamval{T+1,i},1) > 0
                comp{i} = comp{i} + gamval{T+1,i}'*slackval{T+1,i};
                min_comp{i} = min(min_comp{i}, min(gamval{T+1,i}.*slackval{T+1,i}));
            end
            dim{i} = dim{i} + size(gamval{T+1,i},1);
        end
        
        
        
        zvec = new_zvec;
        K = c_K;
        k = c_k;

        
        if params.debug_plot
            plot_fn(xval,z0);
        end
        
        %% update solver stats
        residuals(iter+1) = current_residual;

        %% check for convergence
        
        all_taus_down = true;
        maxtau = 0;
        for i = 1:N
            all_taus_down = all_taus_down && tauval{i}-params.tauval_tolerance*params.tauval_decrease <= params.tauval_tolerance;
            maxtau = max(maxtau, tauval{i});
        end
        
        constraints_satisfied = only_constraint_residual < params.constraint_tolerance;
        
        if iter >= params.max_feasible_iters && constraints_satisfied
            if ((current_residual > params.resid_tolerance && all_taus_down) || (params.open_loop && ~ol_converged))
                disp('Breaking early due to constraint satisfaction and difficulty of progressing');
            end
            if params.debug_print >= 1
                full_optimality = current_residual
                state_optimality = norm(full(state_opt_conds))
                all_taus_down
                ol_converged
                only_constraint_residual
            end
            break;
        elseif (current_residual <= params.resid_tolerance && all_taus_down || current_residual < maxtau)
            if ~params.open_loop && ((~ol_converged && ~fresh_policies) || fresh_policy_count < params.fresh_policy_threshold)
                ol_converged = true;
                if params.debug_print >= 1
                    disp('converged to OL solution');
                end
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
                        full_ineq_jac{t,i} = eval_full_jac_ineq_constraint{t,i}(xu);
                        Sig = sparse(diag(gamval{t,i}./slackval{t,i}));
                        ineq_const_val{t,i} = eval_ineq_constraint{t,i}([xval{t};uval{t,i}]);
                        grad = full(eval_grad_cost{t,i}(xu) +...
                                    full_ineq_jac{t,i}'*(Sig*ineq_const_val{t,i} -...
                                                    gamval{t,i} -...
                                                    tauval{i} * slackval{t,i}.^(-1)));
                        hess = eval_hess_cost{t,i}(xu)+...
                               eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                               eval_hess_dyn_mult{t,i}(xu,lamval{t,i})-...
                               eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t,i});
                        aug = full_ineq_jac{t,i}'*Sig*full_ineq_jac{t,i};
                        hess = full(hess+aug);

                        Q{t,i} = [0 grad'; 
                                  grad hess];
                    end
                end
                for i = 1:N
                    H{T+1,i} = [full(eval_constraint{T+1,i}(xval{T+1})), full(eval_jac_constraint{T+1,i}(xval{T+1}))];

                    full_ineq_jac{T+1,i} = eval_jac_ineq_constraint{T+1,i}(xval{T+1});
                    Sig = sparse(diag(gamval{T+1,i}./slackval{T+1,i}));
                    ineq_const_val{T+1,i} = eval_ineq_constraint{T+1,i}(xval{T+1});
                    grad = full(eval_grad_cost{T+1,i}(xval{T+1}) +...
                                    full_ineq_jac{T+1,i}'*(Sig*ineq_const_val{T+1,i} -...
                                                    gamval{T+1,i} -...
                                                    tauval{i} * slackval{T+1,i}.^(-1)));


                    hess = eval_hess_cost{T+1,i}(xval{T+1})+...
                        eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i})-...
                        eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1,i});
                    aug = full_ineq_jac{T+1,i}'*Sig*full_ineq_jac{T+1,i};
                    hess = full(hess+aug);
                    Q{T+1,i} = [0 grad'; 
                                grad hess];
                end
                [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
                for t = 1:T
                    for i = 1:N
                        if slackval{t,i} < 1e-3
                            for tt = t:-1:max(1,t-3)
                                for j = 1:N
                                    K{tt,j} = zeros(size(K{tt,j}));
                                    k{tt,j} = zeros(size(k{tt,j}));
                                end
                            end
                        end
                    end
                end
                fresh_policies = true
                fresh_policy_count = fresh_policy_count + 1;
                reg = params.dev_penalty * params.dev_decay^(fresh_policy_count-1);
%                 plot_fn(xval,z0);
%                 params.debug_plot = false;
                for i = 1:N
                    tauval{i} = params.tauval_init{i};
                end
%                 restoration = false;
                for t = 1:T
                    for i = 1:N
                        devopt{t,i} = uval{t,i}-(K{t,i}*xval{t}+k{t,i});
                        devval{t,i} = devopt{t,i};
                        % Need to think about what gradient of lagrangian
                        % is for devvals 
                        psival{t,i}(start_m{i}+1:start_m{i}+m{i}) = -(devval{t,i}-devopt{t,i})*params.dev_penalty;
%                         devval{t,i} = zeros(m{i},1);
                    end
                end
                dynamic_conds = eval_dynamics(zvec);
                all_constraints = [dynamic_conds];
                new_necessary_conditions = [dynamic_conds];
                all_conds = [dynamic_conds];
                state_opt_conds = [];
                for i = 1:N
                    % necessary conditions for state and other_control
                    % optimalities need to be corrected to account for
                    % anticipated-intent constraints/multipliers.
                    constraint_opt = eval_constraint_violation{i}(zvec);
                    ineq_constraint_opt = eval_ineq_constraint_violation{i}(zvec);
                    state_opt = eval_state_optimality{i}(zvec);
                    control_opt = eval_control_optimality{i}(zvec);
                    all_control_opt = eval_full_control_optimality{i}(zvec);
                    other_control_opt = eval_other_control_optimality{i}(zvec);
                    complementarity_opt = [];
                    policy_opt = [];
                    deviation_opt = [];

                    ind = 0;
                    ind2 = 0;

                    for t = 1:T
                        c_KK = [];
                        c_kk = [];
                        for j = 1:N
                            c_KK = [c_KK; K{t,j}];
                            c_kk = [c_kk; k{t,j}];
                        end
                        deviation_opt = [deviation_opt; 
                                         reg*(devval{t,i}-devopt{t,i}) + psival{t,i}(start_m{i}+1:start_m{i}+m{i})]; 
                        if t > 1
                            state_opt(ind+1:ind+n) = state_opt(ind+1:ind+n) + c_KK'*psival{t,i};
                            ind = ind+n;
                        end

                        policy_opt = [policy_opt; K{t,i}*xval{t}+k{t,i}+devval{t,i}-uval{t,i}];

                        all_control_opt(ind2+1:ind2+all_m) = all_control_opt(ind2+1:ind2+all_m) - psival{t,i};
                        ind2 = ind2+all_m;
                    end
                    for t= 1:T+1
                        complementarity_opt = [complementarity_opt; gamval{t,i}.*slackval{t,i}-tauval{i}];
                    end
                    private_conds{i} = [constraint_opt;
                                        state_opt;
                                        ineq_constraint_opt;
                                        all_control_opt;
                                        complementarity_opt;
                                        policy_opt;
                                        deviation_opt];
                                    
                    private_fine_resids{i,1} = norm(full(constraint_opt));
                    private_fine_resids{i,2} = norm(full(state_opt));
                    private_fine_resids{i,3} = norm(full(ineq_constraint_opt));
                    private_fine_resids{i,4} = norm(full(raw_ineq_constraint_opt));
                    private_fine_resids{i,5} = norm(full(all_control_opt));
                    private_fine_resids{i,6} = norm(full(complementarity_opt));
                    private_fine_resids{i,7} = norm(full(deviation_opt));
                    private_fine_resids{i,8} = norm(full(policy_opt));     
                    
                    all_constraints = [all_constraints;
                                       min(full(eval_raw_ineq_constraint_violation{i}(new_zvec)),0);
                                       constraint_opt];

                    new_necessary_conditions = [new_necessary_conditions; 
                                                private_conds{i}];
                    all_conds = [all_conds; 
                                state_opt;
                                private_conds{i}];

                    private_conds{i} = [private_conds{i}; 
                                        dynamic_conds];
                    state_opt_conds = [state_opt_conds;
                                        state_opt];
                    new_private_residual{i} = norm(full(private_conds{i}));
                end

                new_current_residual = norm(full(new_necessary_conditions));
                only_constraint_residual = norm(full(all_constraints));
                current_residual = new_current_residual;
                private_residual = new_private_residual;
                
                
%                 restoration = true;
                continue;
            elseif all_taus_down
                state_optimality=norm(full(state_opt_conds))
                break;
            else
                if params.independent_updates
                    for i = 1:N
                        if params.basic_tauschedule
                            if private_residual{i} < tauval{i}+params.tauval_tolerance*params.tauval_decrease && tauval{i}-params.tauval_tolerance*params.tauval_decrease > params.tauval_tolerance
                                tauval{i} = params.tauval_decrease*tauval{i};
                                tauval_decreased = true;
                            end
                        else
                            if ~restoration
                                eta = min_comp{i}/(comp{i}/dim{i});
                                sigma = 0.1*min(0.05*(1-eta)/eta,2)^3;
                                tauval{i} = sigma*comp{i}/dim{i};
                            end
                        end
                    end
                else
                    if ~params.wait_for_policy_to_update_tau || fresh_policy_count >= params.fresh_policy_threshold
                        if params.basic_tauschedule 
                            if(current_residual < tauval{i}+params.tauval_tolerance*params.tauval_decrease)
                                for i = 1:N
                                    tauval{i} = params.tauval_decrease*tauval{i};
                                    tauval_decreased = true;
                                end
                            end
                        else
                            max_tau = -inf;
                            for i = 1:N
                                if dim{i} > 0
                                    eta = min_comp{i}/(comp{i}/dim{i});
                                    sigma = 0.1*min(0.05*(1-eta)/eta,2)^3;
                                    tauval{i} = sigma*comp{i}/dim{i};
                                    max_tau = max(max_tau, tauval{i});
                                end
                            end
                            for i = 1:N
                                tauval{i} = max_tau;
                            end
                        end
                    end
                end
            end
        else
            tauval_decreased = false;
            fresh_policies = false;
            if params.independent_updates
                for i = 1:N
                    if params.basic_tauschedule
                        if private_residual{i} < tauval{i}+params.tauval_tolerance*params.tauval_decrease && tauval{i}-params.tauval_tolerance*params.tauval_decrease > params.tauval_tolerance
                            tauval{i} = params.tauval_decrease*tauval{i};
                            tauval_decreased = true;
                        end
                    else
                        if ~restoration
                            eta = min_comp{i}/(comp{i}/dim{i});
                            sigma = 0.1*min(0.05*(1-eta)/eta,2)^3;
                            tauval{i} = sigma*comp{i}/dim{i};
                        end
                    end
                end
            else
                if ~params.wait_for_policy_to_update_tau || fresh_policy_count >= params.fresh_policy_threshold
                    if params.basic_tauschedule 
                        if(current_residual < tauval{i}+params.tauval_tolerance*params.tauval_decrease)
                            for i = 1:N
                                tauval{i} = params.tauval_decrease*tauval{i};
                                tauval_decreased = true;
                            end
                        end
                    else
                        max_tau = -inf;
                        for i = 1:N
                            if dim{i} > 0
                                eta = min_comp{i}/(comp{i}/dim{i});
                                sigma = 0.1*min(0.05*(1-eta)/eta,2)^3;
                                tauval{i} = sigma*comp{i}/dim{i};
                                max_tau = max(max_tau, tauval{i});
                            end
                        end
                        for i = 1:N
                            tauval{i} = max_tau;
                        end
                    end
                end
            end
%             if tauval_decreased
            tic;
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
                    full_ineq_jac{t,i} = eval_full_jac_ineq_constraint{t,i}(xu);
                    Sig = sparse(diag(gamval{t,i}./slackval{t,i}));
                    ineq_const_val{t,i} = eval_ineq_constraint{t,i}([xval{t};uval{t,i}]);
                    grad = full(eval_grad_cost{t,i}(xu) +...
                                full_ineq_jac{t,i}'*(Sig*ineq_const_val{t,i} -...
                                                gamval{t,i} -...
                                                tauval{i} * slackval{t,i}.^(-1)));
                    hess = eval_hess_cost{t,i}(xu)+...
                           eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                           eval_hess_dyn_mult{t,i}(xu,lamval{t,i})-...
                           eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t,i});
                    aug = full_ineq_jac{t,i}'*Sig*full_ineq_jac{t,i};
                    hess = full(hess+aug);

                    Q{t,i} = [0 grad'; 
                              grad hess];
                end
            end
            for i = 1:N
                H{T+1,i} = [full(eval_constraint{T+1,i}(xval{T+1})), full(eval_jac_constraint{T+1,i}(xval{T+1}))];

                full_ineq_jac{T+1,i} = eval_jac_ineq_constraint{T+1,i}(xval{T+1});
                Sig = sparse(diag(gamval{T+1,i}./slackval{T+1,i}));
                ineq_const_val{T+1,i} = eval_ineq_constraint{T+1,i}(xval{T+1});
                grad = full(eval_grad_cost{T+1,i}(xval{T+1}) +...
                                full_ineq_jac{T+1,i}'*(Sig*ineq_const_val{T+1,i} -...
                                                gamval{T+1,i} -...
                                                tauval{i} * slackval{T+1,i}.^(-1)));


                hess = eval_hess_cost{T+1,i}(xval{T+1})+...
                    eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i})-...
                    eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1,i});
                aug = full_ineq_jac{T+1,i}'*Sig*full_ineq_jac{T+1,i};
                hess = full(hess+aug);
                Q{T+1,i} = [0 grad'; 
                            grad hess];
            end
            if ~params.open_loop && (~params.converge_on_every_policy || ol_converged)
                [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
                if params.converge_on_every_policy
                    ol_converged = false;
                end
                fresh_policies = true;
                fresh_policy_count = fresh_policy_count + 1;
            end
%             end
        end
        deriv_terms_time = toc;
        times = [solve_time,extract_slack_time,ls_time,deriv_terms_time];
        if params.debug_print >= 2
            tauval
            private_residual
            times
        end
    end
    total_iters = iter
    full_sol{1} = xval;
    full_sol{2} = uval;
    full_sol{3} = lamval;
    full_sol{4} = muval;
    full_sol{5} = psival;

end