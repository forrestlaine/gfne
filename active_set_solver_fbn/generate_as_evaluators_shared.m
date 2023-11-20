function evaluators = generate_as_evaluators_shared(f,... % dynamics function
                                             h,... % time-indexed dict of constraint functions
                                             hi,...% time-indexed dict of constraint ownership
                                             g,... % time-indexed dict of inequality constraint functions
                                             gi,...% time-indexed dict of inequality constraint ownership
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
        lam_vars{i} = [];
    end
    mu_vars = [];
    gam_vars = [];
    
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
        end
        constraint{t} = h{t}(states_controls);
        
        gam{t} = MX.sym(['gam_' num2str(t)], numel(g{t}));
        for j = 1:numel(g{t})
            ineq_constraint{t}{j} = g{t}{j}(states_controls);
            jac_ineq_constraint{t}{j} = jacobian(ineq_constraint{t}{j},states_controls);
            evaluators.eval_ineq_constraint{t}{j} = Function('ineq_constraint', {states_controls},{ineq_constraint{t}{j}});
            evaluators.eval_ineq_jac_constraint{t}{j} = Function('jac_ineq_constraint', {states_controls},{jac_ineq_constraint{t}{j}});
        end
        mu{t} = MX.sym(['mu_' num2str(t)], size(constraint{t},1));
        jac_constraint{t} = jacobian(constraint{t},states_controls);
        
        evaluators.eval_constraint{t} = Function('constraint', {states_controls},{constraint{t}});
        evaluators.eval_jac_constraint{t} = Function('jac_constraint', {states_controls},{jac_constraint{t}});
                
        for i = 1:N
            all_vars = [all_vars; u{t,i}; lam{t,i}];
            lam_vars{i} = [lam_vars{i}; lam{t,i}];
            control_vars{i} = [control_vars{i}; u{t,i}];
            all_control_vars = [all_control_vars; u{t,i}];
        end
        all_vars = [all_vars; mu{t}; gam{t}];
        mu_vars = [mu_vars; mu{t}];
        gam_vars = [gam_vars; gam{t}];
        
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
            
            mu_mask = zeros(size(mu{t}));
            gam_mask = zeros(size(gam{t}));
            for j = 1:size(mu{t},1)
                if ismember(i,hi{t}{j})
                    mu_mask(j) = 1;
                end
            end
            for j = 1:size(gam{t},1)
                if ismember(i,gi{t}{j})
                    gam_mask(j) = 1;
                end
            end

            hess_constraint_mult{t,i} = jacobian(jac_constraint{t}'*(mu{t}.*mu_mask),states_controls);
%             hess_ineq_constraint_mult{t,i} = jacobian(jac_ineq_constraint{t}'*(gam{t}.*gam_mask),states_controls);
            hess_dyn_mult{t,i} = jacobian(lam{t,i}'*jac_pred{t},states_controls);
            
            evaluators.eval_cost{t,i} = Function('cost', {states_controls},{cost{t,i}});
            evaluators.eval_grad_cost{t,i} = Function('cost', {states_controls},{grad_cost{t,i}});
            evaluators.eval_hess_cost{t,i} = Function('cost', {states_controls},{hess_cost{t,i}});
            
            evaluators.eval_hess_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,mu{t}},{hess_constraint_mult{t,i}});
%             evaluators.eval_hess_ineq_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,gam{t}},{hess_ineq_constraint_mult{t,i}});
            evaluators.eval_hess_dyn_mult{t,i} = Function('hess_dyn', {states_controls,lam{t,i}},{hess_dyn_mult{t,i}});
          
            
            % DOES NOT INCLUDE INEQUALITIES RIGHT NOW!
            lagrangian{i} = lagrangian{i} + cost{t,i} +...
                (mu{t}.*mu_mask)'*constraint{t} +...
                lam{t,i}'*(pred{t}-x{t+1});% +...
%                 (gam{t}.*gam_mask)'*ineq_constraint{t};
        end

    end
    final_state = x{T+1};
    constraint{T+1} = h{T+1}(final_state);
    
    
    gam{T+1} = MX.sym(['gam_' num2str(T+1)], numel(g{T+1}));
    evaluators.eval_ineq_constraint{T+1} = {};
    evaluators.eval_jac_ineq_constraint{T+1} = {};
    
    for j = 1:numel(g{T+1})
        ineq_constraint{T+1}{j} = g{T+1}{j}(final_state);
        jac_ineq_constraint{T+1}{j} = jacobian(ineq_constraint{T+1}{j},final_state);
        evaluators.eval_ineq_constraint{T+1}{j} = Function('constraint', {final_state},{ineq_constraint{T+1}{j}});
        evaluators.eval_jac_ineq_constraint{T+1}{j} = Function('jac_ineq_constraint', {final_state},{jac_ineq_constraint{T+1}{j}});
    end
    jac_constraint{T+1} = jacobian(constraint{T+1},final_state);
  
    evaluators.eval_constraint{T+1} = Function('constraint', {final_state},{constraint{T+1}});
    evaluators.eval_jac_constraint{T+1} = Function('jac_constraint', {final_state},{jac_constraint{T+1}});

    mu{T+1} = MX.sym(['mu_' num2str(T+1)], size(constraint{T+1},1));
    mu_vars = [mu_vars; mu{T+1}];
    gam_vars = [gam_vars; gam{T+1}];
    all_vars = [all_vars; mu{T+1}; gam{T+1}];
    for i = 1:N

        cost{T+1,i} = l{T+1,i}(final_state);
        
        mu_mask = zeros(size(mu{T+1}));
        gam_mask = zeros(size(gam{T+1}));
        for j = 1:size(mu{T+1},1)
            if ismember(i,hi{T+1}{j})
                mu_mask(j) = 1;
            end
        end
        for j = 1:size(gam{T+1},1)
            if ismember(i,gi{T+1}{j})
                gam_mask(j) = 1;
            end
        end
        
        grad_cost{T+1,i} = gradient(cost{T+1,i},final_state);
        hess_cost{T+1,i} = jacobian(grad_cost{T+1,i},final_state);
        
        hess_constraint_mult{T+1,i} = jacobian(jac_constraint{T+1}'*(mu{T+1}.*mu_mask),final_state);
%         hess_ineq_constraint_mult{T+1,i} = jacobian(jac_ineq_constraint{T+1}'*(gam{T+1}.*gam_mask),final_state);
        
        
        evaluators.eval_cost{T+1,i} = Function('cost', {final_state},{cost{T+1,i}});
        evaluators.eval_grad_cost{T+1,i} = Function('cost', {final_state},{grad_cost{T+1,i}});
        evaluators.eval_hess_cost{T+1,i} = Function('cost', {final_state},{hess_cost{T+1,i}});
        evaluators.eval_hess_constraint_mult{T+1,i} = Function('hess_constraint_mult', {final_state,mu{T+1}},{hess_constraint_mult{T+1,i}});
%         evaluators.eval_hess_ineq_constraint_mult{T+1,i} = Function('hess_ineq_constraint_mult', {final_state,gam{T+1}},{hess_ineq_constraint_mult{T+1,i}});
        
        % DOES NOT INCLUDE INEQUALITIES
        lagrangian{i} = lagrangian{i} + cost{T+1,i} + (mu{T+1}.*mu_mask)'*constraint{T+1};% - (gam{T+1}.*gam_mask)'*ineq_constraint{T+1};
    end
    
    % NOT IMPLEMENTED: HELPER EVALUATORS FOR MERIT FUNCTION
    
    
%     evaluators.eval_dynamics = Function('dyn', {all_vars},{shared_constraints});
%     for i = 1:N
%         constraint_violation{i} = gradient(lagrangian{i},mu_vars{i}); % h(x)
%         raw_ineq_constraint_violation{i} = -gradient(lagrangian{i},gam_vars{i});
%         
%         evaluators.eval_constraint_violation{i} = Function('con', {all_vars}, {constraint_violation{i}});
%         evaluators.eval_raw_ineq_constraint_violation{i} = Function('con', {all_vars}, {raw_ineq_constraint_violation{i}});
%         state_optimality{i} = gradient(lagrangian{i},state_vars);
%         evaluators.eval_state_optimality{i} = Function('con', {all_vars}, {state_optimality{i}});
%         control_optimality{i} = gradient(lagrangian{i},control_vars{i});
%         full_control_optimality{i} = gradient(lagrangian{i},all_control_vars);
%         evaluators.eval_control_optimality{i} = Function('con', {all_vars}, {control_optimality{i}});
%         evaluators.eval_full_control_optimality{i} = Function('con', {all_vars}, {full_control_optimality{i}});
%         other_control_optimality{i} = gradient(lagrangian{i},other_control_vars{i});
%         evaluators.eval_other_control_optimality{i} = Function('con', {all_vars}, {other_control_optimality{i}});
%     end
    
end