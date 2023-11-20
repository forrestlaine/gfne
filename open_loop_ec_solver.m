%% Iterative EC solver

function [residuals, xval] = open_loop_ec_solver(f,... % dynamics function
                         h,... % time-indexed dict of player-indexed equality constraint functions
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
        gam_vars{i} = [];
        lam_vars{i} = [];
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
            constraint{t,i} = h{t,i}([x{t};u{t,i}]);
            iconstraint{t,i} = g{t,i}([x{t};u{t,i}]);
            mu{t,i} = MX.sym(['mu_' num2str(t) '_' num2str(i)], size(constraint{t,i},1));
            gam{t,i} = MX.sym(['gam_' num2str(t) '_' num2str(i)], size(iconstraint{t,i},1));
            slack{t,i} = MX.sym(['slack_' num2str(t) '_' num2str(i)], size(iconstraint{t,i},1));
            
            jac_constraint{t,i} = jacobian(constraint{t,i},[x{t};u{t,i}]);
            jac_iconstraint{t,i} = jacobian(iconstraint{t,i},[x{t};u{t,i}]);
                        
            eval_constraint{t,i} = Function('constraint', {[x{t};u{t,i}]},{constraint{t,i}});
            eval_iconstraint{t,i} = Function('iconstraint', {[x{t};u{t,i}]},{iconstraint{t,i}});
            eval_jac_constraint{t,i} = Function('jac_constraint', {[x{t};u{t,i}]},{jac_constraint{t,i}});
            eval_jac_iconstraint{t,i} = Function('jac_iconstraint', {[x{t};u{t,i}]},{jac_iconstraint{t,i}});
            
            all_vars = [all_vars; u{t,i}; mu{t,i}; gam{t,i}; lam{t,i}];
            lam_vars{i} = [lam_vars{i}; lam{t,i}];
            control_vars{i} = [control_vars{i}; u{t,i}];
            mu_vars{i} = [mu_vars{i}; mu{t,i}];
            gam_vars{i} = [gam_vars{i}; gam{t,i}];
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
            full_jac_iconstraint{t,i} = jacobian(iconstraint{t,i}, states_controls);
            hess_constraint_mult{t,i} = jacobian(full_jac_constraint{t,i}'*mu{t,i},states_controls);
            hess_iconstraint_mult{t,i} = jacobian(full_jac_iconstraint{t,i}'*gam{t,i},states_controls);
            eval_hess_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,mu{t,i}},{hess_constraint_mult{t,i}});
            eval_hess_iconstraint_mult{t,i} = Function('hess_iconstraint_mult', {states_controls,gam{t,i}},{hess_iconstraint_mult{t,i}});
            
            hess_dyn_mult{t,i} = jacobian(lam{t,i}'*jac_pred{t},states_controls);
            
            eval_cost{t,i} = Function('cost', {states_controls},{cost{t,i}});
            eval_grad_cost{t,i} = Function('cost', {states_controls},{grad_cost{t,i}});
            eval_hess_cost{t,i} = Function('cost', {states_controls},{hess_cost{t,i}});
            
            eval_hess_dyn_mult{t,i} = Function('hess_dyn', {states_controls,lam{t,i}},{hess_dyn_mult{t,i}});
            lag{t,i} = cost{t,i}+lam{t,i}'*pred{t};
            jac_lag{t,i} = gradient(lag{t,i}, states_controls);
            hess_lag{t,i} = jacobian(jac_lag{t,i}, states_controls);
            eval_hess_lag{t,i} = Function('laghess', {states_controls,lam{t,i}},{hess_lag{t,i}});
            
            lagrangian{i} = lagrangian{i} + cost{t,i} + mu{t,i}'*constraint{t,i} + gam{t,i}'*iconstraint{t,i} + lam{t,i}'*(pred{t}-x{t+1});
        end

    end
    final_state = x{T+1};
    for i = 1:N
        constraint{T+1,i} = h{T+1,i}(final_state);
        iconstraint{T+1,i} = g{T+1,i}(final_state);
        mu{T+1,i} = MX.sym(['mu_' num2str(T+1) '_' num2str(i)], size(constraint{T+1,i},1));
        gam{T+1,i} = MX.sym(['gam_' num2str(T+1) '_' num2str(i)], size(iconstraint{T+1,i},1));
        slack{T+1,i} = MX.sym(['slack_' num2str(T+1) '_' num2str(i)], size(iconstraint{T+1,i},1));
        
        all_vars = [all_vars; mu{T+1,i}; gam{T+1,i}];
        mu_vars{i} = [mu_vars{i}; mu{T+1,i}];
        gam_vars{i} = [gam_vars{i}; gam{T+1,i}];
        
        cost{T+1,i} = l{T+1,i}(final_state);
        grad_cost{T+1,i} = gradient(cost{T+1,i},final_state);
        hess_cost{T+1,i} = jacobian(grad_cost{T+1,i},final_state);
        
        jac_constraint{T+1,i} = jacobian(constraint{T+1,i},final_state);
        hess_constraint_mult{T+1,i} = jacobian(jac_constraint{T+1,i}'*mu{T+1,i},final_state);
        eval_constraint{T+1,i} = Function('constraint', {final_state},{constraint{T+1,i}});
        eval_jac_constraint{T+1,i} = Function('jac_constraint', {final_state},{jac_constraint{T+1,i}});
        eval_hess_constraint_mult{T+1,i} = Function('hess_constraint_mult', {final_state,mu{T+1,i}},{hess_constraint_mult{T+1,i}});
        
        jac_iconstraint{T+1,i} = jacobian(iconstraint{T+1,i},final_state);
        hess_iconstraint_mult{T+1,i} = jacobian(jac_iconstraint{T+1,i}'*gam{T+1,i},final_state);
        eval_iconstraint{T+1,i} = Function('iconstraint', {final_state},{iconstraint{T+1,i}});
        eval_jac_iconstraint{T+1,i} = Function('jac_iconstraint', {final_state},{jac_iconstraint{T+1,i}});
        eval_hess_iconstraint_mult{T+1,i} = Function('hess_iconstraint_mult', {final_state,gam{T+1,i}},{hess_iconstraint_mult{T+1,i}});
        
        eval_cost{T+1,i} = Function('cost', {final_state},{cost{T+1,i}});
        eval_grad_cost{T+1,i} = Function('cost', {final_state},{grad_cost{T+1,i}});
        eval_hess_cost{T+1,i} = Function('cost', {final_state},{hess_cost{T+1,i}});
        
        lagrangian{i} = lagrangian{i} + cost{T+1,i} + mu{T+1,i}'*constraint{T+1,i} + gam{T+1,i}'*iconstraint{T+1,i};
    end
    
    necessary_conditions = [shared_constraints];
    eval_dynamics = Function('dyn', {all_vars},{necessary_conditions});
    full_vars = [state_vars];
    for i = 1:N
        constraint_violation{i} = gradient(lagrangian{i},mu_vars{i});
        eval_constraint_violation{i} = Function('con', {all_vars}, {constraint_violation{i}});
        state_optimality{i} = gradient(lagrangian{i},state_vars);
        eval_state_optimality{i} = Function('con', {all_vars}, {state_optimality{i}});
        control_optimality{i} = gradient(lagrangian{i},control_vars{i});
        eval_control_optimality{i} = Function('con', {all_vars}, {control_optimality{i}});
        other_control_optimality{i} = gradient(lagrangian{i},other_control_vars{i});
        eval_other_control_optimality{i} = Function('con', {all_vars}, {other_control_optimality{i}});
        necessary_conditions = [necessary_conditions; constraint_violation{i}; state_optimality{i}; control_optimality{i}];
        full_vars = [full_vars; control_vars{i}; lam_vars{i}; mu_vars{i}];
    end
    lag_hess = jacobian(necessary_conditions, full_vars);
    eval_lag_hess = Function('laghess', {all_vars}, {lag_hess});
    eval_necessary_conditions = Function('FONCs', {all_vars}, {necessary_conditions}); 
    %% initialize solution from z0 using no control
    xval{1} = z0;
    for t = 1:T
        xu = [xval{t}];
        for i = 1:N
            uval{t,i} = zeros(m{i},1);
            lamval{t,i} = zeros(n,1);
            muval{t,i} = zeros(size(h{t,i}([xval{t};uval{t,i}]),1),1);
            if t<T
                psival{t,i} = zeros(all_m-m{i},1);
            end
            xu = [xu; uval{t,i}];
        end
        xval{t+1} = f(xu);
    end
    
    for i = 1:N
        muval{T+1,i} = zeros(size(h{T+1,i}(xval{T+1}),1),1);
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
            hess = eval_hess_cost{t,i}(xu)+...
                   eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                   eval_hess_dyn_mult{t,i}(xu,lamval{t,i});
            Q{t,i} = [0 full(eval_grad_cost{t,i}(xu))'; 
                full(eval_grad_cost{t,i}(xu)) full(hess)];
        end
    end
    for i = 1:N
        H{T+1,i} = [full(eval_constraint{T+1,i}(xval{T+1})), full(eval_jac_constraint{T+1,i}(xval{T+1}))];
        hess = eval_hess_cost{T+1,i}(xval{T+1})+ eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i});
        Q{T+1,i} = [0 full(eval_grad_cost{T+1,i}(xval{T+1}))'; 
                      full(eval_grad_cost{T+1,i}(xval{T+1})) full(hess)];
    end
%     %% initialize policies

    for t = 1:T
        for i = 1:N
            K{t,i} = zeros(m{i},n);
        end
    end
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
    end
    current_residual = norm(full(cur_necessary_conditions));
    residuals(1) = current_residual; 
    prev_alpha = 1;
    prev_s = 1;
    %% iteratively solve for solution
    for iter = 1:100 % max iterations
        %% extract update / solve for search direction
        [dX,dU,dL,dM,dP] = solve_ec_lq_game_d(F,H,Q,N,T,K,zeros(size(z0)));
        
            
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
            for i = 1:N
                % necessary conditions for state and other_control
                % optimalities need to be corrected to account for
                % anticipated-intent constraints/multipliers.
                constraint_opt = eval_constraint_violation{i}(new_zvec);
                state_opt = eval_state_optimality{i}(new_zvec);
                control_opt = eval_control_optimality{i}(new_zvec);

  
                
                if s == num_linesearch_iters
                    disp('check');
                end

                new_necessary_conditions = [new_necessary_conditions; 
                                            state_opt;
                                            constraint_opt; 
                                            control_opt];
            end

            new_current_residual = norm(full(new_necessary_conditions));
            ls_residual = new_current_residual;
            ls_resids(s+1) = new_current_residual;

            if (new_current_residual < current_residual)
                disp('linesearch accepted!');
                current_residual = new_current_residual;
                prev_alpha = alpha/beta;
                prev_s = s-1;
                linesearch_succeeded = true;
                break;
            else
                alpha = alpha * beta;
            end
        end

        %% Accept values after step
        if linesearch_succeeded
            H = c_H;
            Q = c_Q;
            F = c_F;
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
        else
            disp('ls did not converge!');
            break;
                           
        end
        
        %% update solver stats
        residuals(iter+1) = current_residual;

        %% check for convergence
        current_residual
        if current_residual < 1e-3
            break;
        end
    end
    
    close all;
    figure; hold on;
    xx(1,:) = z0;
    spot_a = plot(xx(1,1), xx(1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    spot_b = plot(xx(1,5), xx(1,6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    pause(0.1);
    for t = 1:T
        xx(t+1,:) = xval{t+1};
        spot_a = plot(xx(t+1,1), xx(t+1,2),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        spot_b = plot(xx(t+1,5), xx(t+1,6),'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        axis([-10 10 -10 10]);
        pause(0.1);
    end
    disp('plotted');
    
    sol = current_residual;
end