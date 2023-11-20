%% Iterative EC solver

function [residuals, xval] = ego_centric_ec_solver(f,... % dynamics function
                         h,... % time-indexed dict of player-indexed constraint functions
                         l,... % time-indexed dict of player-indexed cost functions
                         n,... % dimension of shared state
                         m,... % player-indexed dict of player input dimensions
                         N,... % number of players
                         T,... % time steps
                         p,... % ego index
                         z0)   % initial state
    open_loop = false;
    beta = 0.5;
    rho = 0.5;
    eta = 0.1;
    failed_linesearches = 0;
    import casadi.*
    %% initialize variables and jacobian/hessian functions
    
    shared_constraints = [];
    x{1} = MX.sym(['x_' num2str(1)], n);
    all_vars = [x{1}];
    state_vars = [x{1}];
    all_m = 0;
    fp = 0;
    for i = 1:N
        lam{1,i} = MX.sym(['lam_' num2str(1) '_' num2str(i)],n);
        lagrangian{i} = lam{1,i}'*(z0-x{1});
        all_m = all_m + m{i};
        control_vars{i} = [];
        other_control_vars{i} = [];
        mu_vars{i} = [];
        lam_vars{i} = [lam{1,i}];
        all_vars = [all_vars; lam{1,i}];
    end
    

    
    for t = 1:T
        x{t+1} = MX.sym(['x_' num2str(t+1)], n);
        states_controls = [x{t}];
        
        all_vars = [all_vars; x{t+1}];
        state_vars = [state_vars; x{t+1}];

        
        % Constraints
        for i = 1:N
            u{t,i} = MX.sym(['u_' num2str(t) '_' num2str(i)], m{i});
            lam{t+1,i} = MX.sym(['lam_' num2str(t) '_' num2str(i)], n);
            states_controls = [states_controls; u{t,i}];
            constraint{t,i} = h{t,i}([x{t};u{t,i}]);
            mu{t,i} = MX.sym(['mu_' num2str(t) '_' num2str(i)], size(constraint{t,i},1));
            jac_constraint{t,i} = jacobian(constraint{t,i},[x{t};u{t,i}]);
                        
            eval_constraint{t,i} = Function('constraint', {[x{t};u{t,i}]},{constraint{t,i}});
            eval_jac_constraint{t,i} = Function('jac_constraint', {[x{t};u{t,i}]},{jac_constraint{t,i}});
            
            all_vars = [all_vars; u{t,i}; mu{t,i}; lam{t+1,i}];
            lam_vars{i} = [lam_vars{i}; lam{t+1,i}];
            control_vars{i} = [control_vars{i}; u{t,i}];
            mu_vars{i} = [mu_vars{i}; mu{t,i}];
        end
        
 
        for i = 1:N
            for j = 1:N
                if j ~= i
                    other_control_vars{i} = [other_control_vars{i}; u{t,j}];
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
            if i == p
                fp = fp + cost{t,p};
            end
            grad_cost{t,i} = gradient(cost{t,i},states_controls);
            hess_cost{t,i} = jacobian(grad_cost{t,i},states_controls);
            
            
            full_jac_constraint{t,i} = jacobian(constraint{t,i}, states_controls);
            hess_constraint_mult{t,i} = jacobian(full_jac_constraint{t,i}'*mu{t,i},states_controls);
            hess_dyn_mult{t,i} = jacobian(lam{t+1,i}'*jac_pred{t},states_controls);
            
            eval_cost{t,i} = Function('cost', {states_controls},{cost{t,i}});
            eval_grad_cost{t,i} = Function('cost', {states_controls},{grad_cost{t,i}});
            eval_hess_cost{t,i} = Function('cost', {states_controls},{hess_cost{t,i}});
            
            eval_hess_constraint_mult{t,i} = Function('hess_constraint_mult', {states_controls,mu{t,i}},{hess_constraint_mult{t,i}});
            eval_hess_dyn_mult{t,i} = Function('hess_dyn', {states_controls,lam{t+1,i}},{hess_dyn_mult{t,i}});
            lag{t,i} = cost{t,i}+lam{t,i}'*pred{t};
            jac_lag{t,i} = gradient(lag{t,i}, states_controls);
            hess_lag{t,i} = jacobian(jac_lag{t,i}, states_controls);
            eval_hess_lag{t,i} = Function('laghess', {states_controls,lam{t,i}},{hess_lag{t,i}});
            
            lagrangian{i} = lagrangian{i} + cost{t,i} + mu{t,i}'*constraint{t,i} + lam{t+1,i}'*(pred{t}-x{t+1});
        end

    end
    final_state = x{T+1};
    for i = 1:N
        constraint{T+1,i} = h{T+1,i}(final_state);
        mu{T+1,i} = MX.sym(['mu_' num2str(T+1) '_' num2str(i)], size(constraint{T+1,i},1));
        all_vars = [all_vars; mu{T+1,i}];
        mu_vars{i} = [mu_vars{i}; mu{T+1,i}];
        cost{T+1,i} = l{T+1,i}(final_state);
        if i == p
            fp = fp + cost{T+1,p};
        end
        jac_constraint{T+1,i} = jacobian(constraint{T+1,i},final_state);
        grad_cost{T+1,i} = gradient(cost{T+1,i},final_state);
        hess_cost{T+1,i} = jacobian(grad_cost{T+1,i},final_state);
        hess_constraint_mult{T+1,i} = jacobian(jac_constraint{T+1,i}'*mu{T+1,i},final_state);
        
        eval_constraint{T+1,i} = Function('constraint', {final_state},{constraint{T+1,i}});
        eval_jac_constraint{T+1,i} = Function('jac_constraint', {final_state},{jac_constraint{T+1,i}});
        eval_cost{T+1,i} = Function('cost', {final_state},{cost{T+1,i}});
        eval_grad_cost{T+1,i} = Function('cost', {final_state},{grad_cost{T+1,i}});
        eval_hess_cost{T+1,i} = Function('cost', {final_state},{hess_cost{T+1,i}});
        eval_hess_constraint_mult{T+1,i} = Function('hess_constraint_mult', {final_state,mu{T+1,i}},{hess_constraint_mult{T+1,i}});
        
        lagrangian{i} = lagrangian{i} + cost{T+1,i} + mu{T+1,i}'*constraint{T+1,i};
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
%     lag_hess = jacobian(necessary_conditions, full_vars);
%     eval_lag_hess = Function('laghess', {all_vars}, {lag_hess});
    eval_necessary_conditions = Function('FONCs', {all_vars}, {necessary_conditions});
    all_primals = [state_vars];
    for i = 1:N
        all_primals = [all_primals; control_vars{i}];
    end
    grad_fp = gradient(fp, all_primals);
    eval_fp = Function('fp', {all_vars}, {fp});
    eval_grad_fp = Function('gfp', {all_vars}, {grad_fp});
    %% initialize solution from z0 using no control
    xval{1} = z0;
    for t = 1:T
        xu = [xval{t}];
        for i = 1:N
            uval{t,i} = zeros(m{i},1);
            lamval{t,i} = zeros(n,1);
            muval{t,i} = zeros(size(h{t,i}([xval{t};uval{t,i}]),1),1);
            psival{t,i} = zeros(all_m-m{i},1);
            xu = [xu; uval{t,i}];
        end
        xval{t+1} = f(xu);
    end
    
    for i = 1:N
        muval{T+1,i} = zeros(size(h{T+1,i}(xval{T+1}),1),1);
        lamval{T+1,i} = zeros(n,1);
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
            hess = eval_hess_cost{t,i}(xu);%+...
%                    eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
%                    eval_hess_dyn_mult{t,i}(xu,lamval{t,i});
            HK{t,i} = 0*eye(n);
            Q{t,i} = [0 full(eval_grad_cost{t,i}(xu))'; 
                full(eval_grad_cost{t,i}(xu)) full(hess)];
        end
    end
    for i = 1:N
        H{T+1,i} = [full(eval_constraint{T+1,i}(xval{T+1})), full(eval_jac_constraint{T+1,i}(xval{T+1}))];
        hess = eval_hess_cost{T+1,i}(xval{T+1});%+ eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i});
        Q{T+1,i} = [0 full(eval_grad_cost{T+1,i}(xval{T+1}))'; 
                      full(eval_grad_cost{T+1,i}(xval{T+1})) full(hess)];
    end
%     %% initialize policies
    [K, k] = solve_ec_lq_game_ol(F,H,Q,N,T);
    for t = 1:T
        for i = 1:N
            KK{t,i} = [];
            kk{t,i} = [];
            for j = 1:N
                if j ~= i
                    KK{t,i} = [KK{t,i}; K{t,j}];
                    kk{t,i} = [kk{t,i}; k{t,j}];
                end
            end
        end
    end
    %% pack initial z0
    zvec = [z0];
    for i = 1:N
        zvec = [zvec; lamval{1,i}];
    end
    for t = 1:T
        zvec = [zvec; xval{t+1}];
        for i = 1:N
            zvec = [zvec; uval{t,i}; muval{t,i}; lamval{t+1,i}];
        end
    end
    for i = 1:N
        zvec = [zvec; muval{T+1,i}];
    end
    current_residual = inf;
    residuals(1) = current_residual; 
    prev_alpha = 1;
    prev_s = 1;
    muk = 1;
    %% iteratively solve for solution
    for iter = 1:100 % max iterations
        %% extract update / solve for search direction
        [dX,dU,dL,dM,dP,sn] = solve_ec_single_game(F,H,Q,N,T,KK,kk,p,zeros(size(z0)));
        step_size = sn
        if step_size < 1e-2
            break;
        end
            
        %% linesearch (update this so we're not packing and upacking data)
        
        alpha = min(1,prev_alpha);   
%         alpha = 1;
        start_s = max(1,prev_s);
%         start_s = 1;
        num_linesearch_iters = 10;

        
        linesearch_succeeded = false;
        
        pk = [];
        all_constraints = [];
        for t = 1:T+1
            pk = [pk;dX{t}];
        end
        for i = 1:N
            for t = 1:T
                pk = [pk; dU{t,i}];
            end
        end
        all_constraints = [z0-xval{1}];
        for t = 1:T
            all_constraints = [all_constraints; F{t}(:,1)];
        end
        for t = 1:T
            all_constraints = [all_constraints; kk{t,p}];
        end
        for t = 1:T+1
            all_constraints = [all_constraints; H{t,p}(:,1)];
        end
        grad_fp = full(eval_grad_fp(zvec)'*pk);
        constraint_norm = norm(full(all_constraints),1)
        muk = max(muk,grad_fp/((1-rho)*constraint_norm));
        directional_derivative = grad_fp - muk*constraint_norm;
        if iter == 1
            current_merit = inf;
        else
            current_merit = full(eval_fp(zvec)) + muk*constraint_norm;
        end
        ls_resids(1) = current_merit;

        for s = start_s:num_linesearch_iters % max linesearch iterations
            % pack candidate values ('c_') and jacobians, hessians, etc.
            c_xval{1} = z0;
            %%
 
            c_lamval{1,p} = lamval{1,p}+alpha*(dL{1,p}-lamval{1,p});

            for t = 1:T
                for i = 1:N
                    c_uval{t,i} = uval{t,i}+alpha*dU{t,i};
                    if i == p
                        c_lamval{t+1,i} = lamval{t+1,i}+alpha*(dL{t+1,i}-lamval{t+1,i});
                        c_muval{t,i} = muval{t,i}+alpha*(dM{t,i}-muval{t,i}); 
                        c_psival{t,i} = psival{t,i}+alpha*(dP{t,i}-psival{t,i});
                    end
                end
                c_xval{t+1} = xval{t+1}+alpha*dX{t+1};

                uu = [];
                for i = 1:N
                    c_H{t,i} = [full(eval_constraint{t,i}([c_xval{t};c_uval{t,i}])),...
                                full(eval_jac_constraint{t,i}([c_xval{t};c_uval{t,i}]))];
                    uu = [uu; c_uval{t,i}];
                end
                xu = [c_xval{t};uu];
                c_F{t} = full([full(eval_pred{t}(xu))-c_xval{t+1}, full(eval_jac_pred{t}(xu))]);

                for i = 1:N
                    if false
                        hess = eval_hess_cost{t,i}(xu)+ ...
                               eval_hess_constraint_mult{t,i}(xu, c_muval{t,i})+ ...
                               eval_hess_dyn_mult{t,i}(xu,c_lamval{t+1,i});
                    else
                        hess = eval_hess_cost{t,i}(xu);
                    end
                    hess(1+1:1+n,1+1:1+n) = hess(1+1:1+n,1+1:1+n);

                    c_Q{t,i} = [0 full(eval_grad_cost{t,i}(xu))'; 
                                full(eval_grad_cost{t,i}(xu)) full(hess)];
                end

            end
            for i = 1:N
                if i == p
                    c_muval{T+1,i} = muval{T+1,i}+alpha*(dM{T+1,i}-muval{T+1,i});
                end
                hess = eval_hess_cost{T+1,i}(c_xval{T+1});% + eval_hess_constraint_mult{T+1,i}(c_xval{T+1},c_muval{T+1,i});
                c_H{T+1,i} = [full(eval_constraint{T+1,i}(c_xval{T+1})), full(eval_jac_constraint{T+1,i}(c_xval{T+1}))];
                c_Q{T+1,i} = [0 full(eval_grad_cost{T+1,i}(c_xval{T+1}))'; 
                            full(eval_grad_cost{T+1,i}(c_xval{T+1})) full(hess)];
            end

            % Note these policy jac terms are evaluated after each candidate step

%             [c_K, c_k,~,~] = solve_ec_lq_game_b(c_F,c_H,c_Q,N,T);
            [c_K, c_k] = solve_ec_lq_game_ol(c_F,c_H,c_Q,N,T);
            for t = 1:T
                for i = 1:N
                    c_KK{t,i} = [];
                    c_kk{t,i} = [];
                    for j = 1:N
                        if j ~= i
                            c_KK{t,i} = [c_KK{t,i}; c_K{t,j}];
                            c_kk{t,i} = [c_kk{t,i}; c_k{t,j}];
                        end
                    end
                end
            end

            new_zvec = [z0];
            for i = 1:N
                if i == p
                    new_zvec = [new_zvec; c_lamval{1,i}];
                else
                    new_zvec = [new_zvec; lamval{1,i}];
                end
            end
            for t = 1:T
                new_zvec = [new_zvec; 
                            c_xval{t+1}];
                for i = 1:N
                    if i == p
                        new_zvec = [new_zvec; 
                                    c_uval{t,i}; 
                                    c_muval{t,i}; 
                                    c_lamval{t+1,i}];
                    else
                        new_zvec = [new_zvec; 
                                    c_uval{t,i}; 
                                    muval{t,i}; 
                                    lamval{t+1,i}];
                    end
                end
            end
            for i = 1:N
                if i == p
                    new_zvec = [new_zvec; 
                                c_muval{T+1,i}];
                else
                    new_zvec = [new_zvec;
                                muval{T+1,i}];
                end
            end


            %% account for feedback constraints in necessary conditions

            all_constraints = [z0-c_xval{1}];
            for t = 1:T
                all_constraints = [all_constraints; c_F{t}(:,1)];
            end
            for t = 1:T
                all_constraints = [all_constraints; c_kk{t,p}];
            end
            for t = 1:T+1
                all_constraints = [all_constraints; c_H{t,p}(:,1)];
            end
            
            if s == num_linesearch_iters
                disp('check');
            end
            
            cost = full(eval_fp(new_zvec));
            constraint_norm = norm(full(all_constraints),1);
            new_merit = cost + muk*constraint_norm;

            ls_resids(s+1) = new_merit;

            if (new_merit < current_merit+alpha*eta*directional_derivative)
                disp('linesearch accepted!');
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
            KK = c_KK;
            kk = c_kk;
            for t = 1:T
                for i = 1:N
                    uval{t,i} = c_uval{t,i};
                end
                lamval{t,p} = c_lamval{t,p};
                muval{t,p} = c_muval{t,p}; 
                psival{t,p} = c_psival{t,p};
                xval{t+1} = c_xval{t+1};
            end
            muval{T+1,p} = c_muval{T+1,p};
            lamval{T+1,p} = c_lamval{T+1,p};

            zvec = new_zvec;
        else
            disp('ls did not converge!');
            break;
        end
        
        %% update solver stats
        residuals(iter+1) = current_merit;
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
        current_merit
        if current_residual < 1e-3
            break;
        end
        if failed_linesearches > 10
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