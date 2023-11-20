function evaluators = generate_evaluators(f,... % dynamics function
                                          h,... % time-indexed dict of player-indexed constraint functions
                                          g,... % time-indexed dict of player-indexed inequality constraint functions
                                          l,... % time-indexed dict of player-indexed cost functions
                                          n,... % dimension of shared state
                                          m,... % player-indexed dict of player input dimensions
                                          N,... % number of players
                                          T)
    import casadi.*                               
                                     
    shared_constraints = [];
    for i = 1:N
        lagrangian{i} = 0;
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
                        
            evaluators.eval_constraint{t,i} = Function('constraint', {state_single_control},{constraint{t,i}});
            evaluators.eval_jac_constraint{t,i} = Function('jac_constraint', {state_single_control},{jac_constraint{t,i}});
            evaluators.eval_ineq_constraint{t,i} = Function('ineq_constraint', {state_single_control},{ineq_constraint{t,i}});
            evaluators.eval_ineq_jac_constraint{t,i} = Function('jac_ineq_constraint', {state_single_control},{jac_ineq_constraint{t,i}});
            
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
        
        evaluators.eval_pred{t} = Function('dynamics', {states_controls}, {pred{t}});
        evaluators.eval_jac_pred{t} = Function('dynamics_jac', {states_controls}, {jac_pred{t}});
        
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
            
            evaluators.eval_cost{t,i} = Function('cost', {states_controls},{cost{t,i}});
            evaluators.eval_grad_cost{t,i} = Function('cost', {states_controls},{grad_cost{t,i}});
            evaluators.eval_hess_cost{t,i} = Function('cost', {states_controls},{hess_cost{t,i}});
            
            evaluators.eval_full_jac_ineq_constraint{t,i} = Function('full_ineq_jac', {states_controls},{full_jac_ineq_constraint{t,i}});
            evaluators.eval_hess_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,mu{t,i}},{hess_constraint_mult{t,i}});
            evaluators.eval_hess_ineq_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,gam{t,i}},{hess_ineq_constraint_mult{t,i}});
            evaluators.eval_hess_dyn_mult{t,i} = Function('hess_dyn', {states_controls,lam{t,i}},{hess_dyn_mult{t,i}});
          
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
        
        evaluators.eval_constraint{T+1,i} = Function('constraint', {final_state},{constraint{T+1,i}});
        evaluators.eval_ineq_constraint{T+1,i} = Function('constraint', {final_state},{ineq_constraint{T+1,i}});
        evaluators.eval_jac_constraint{T+1,i} = Function('jac_constraint', {final_state},{jac_constraint{T+1,i}});
        evaluators.eval_jac_ineq_constraint{T+1,i} = Function('jac_ineq_constraint', {final_state},{jac_ineq_constraint{T+1,i}});
        evaluators.eval_cost{T+1,i} = Function('cost', {final_state},{cost{T+1,i}});
        evaluators.eval_grad_cost{T+1,i} = Function('cost', {final_state},{grad_cost{T+1,i}});
        evaluators.eval_hess_cost{T+1,i} = Function('cost', {final_state},{hess_cost{T+1,i}});
        evaluators.eval_hess_constraint_mult{T+1,i} = Function('hess_constraint_mult', {final_state,mu{T+1,i}},{hess_constraint_mult{T+1,i}});
        evaluators.eval_hess_ineq_constraint_mult{T+1,i} = Function('hess_ineq_constraint_mult', {final_state,gam{T+1,i}},{hess_ineq_constraint_mult{T+1,i}});
        
        lagrangian{i} = lagrangian{i} + cost{T+1,i} + mu{T+1,i}'*constraint{T+1,i} - gam{T+1,i}'*ineq_constraint{T+1,i};
    end
    

    evaluators.eval_dynamics = Function('dyn', {all_vars},{shared_constraints});
    for i = 1:N
        constraint_violation{i} = gradient(lagrangian{i},mu_vars{i}); % h(x)
        raw_ineq_constraint_violation{i} = -gradient(lagrangian{i},gam_vars{i});
        ineq_constraint_violation{i} = -slack_vars{i}-gradient(lagrangian{i},gam_vars{i}); % g(x) - s
        
        evaluators.eval_constraint_violation{i} = Function('con', {all_vars}, {constraint_violation{i}});
        evaluators.eval_ineq_constraint_violation{i} = Function('con', {all_vars}, {ineq_constraint_violation{i}});
        evaluators.eval_raw_ineq_constraint_violation{i} = Function('con', {all_vars}, {raw_ineq_constraint_violation{i}});
        state_optimality{i} = gradient(lagrangian{i},state_vars);
        evaluators.eval_state_optimality{i} = Function('con', {all_vars}, {state_optimality{i}});
        control_optimality{i} = gradient(lagrangian{i},control_vars{i});
        full_control_optimality{i} = gradient(lagrangian{i},all_control_vars);
        evaluators.eval_control_optimality{i} = Function('con', {all_vars}, {control_optimality{i}});
        evaluators.eval_full_control_optimality{i} = Function('con', {all_vars}, {full_control_optimality{i}});
        other_control_optimality{i} = gradient(lagrangian{i},other_control_vars{i});
        evaluators.eval_other_control_optimality{i} = Function('con', {all_vars}, {other_control_optimality{i}});
    end
    
end