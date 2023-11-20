%% Iterative EC solver

function [x_opt,x,duration,eps_opt] = ai_solver(f,... % dynamics function
                                       h,... % time-indexed dict of player-indexed equality constraint functions
                                       g,... % time-indexed dict of player-indexed inequality constraint functions
                                       l,... % time-indexed dict of player-indexed cost functions
                                       n,... % dimension of shared state
                                       m,... % player-indexed dict of player input dimensions
                                       N,... % number of players
                                       T,... % time steps
                                       p,... % ego player index
                                       z0,...
                                       varargin)   % initial state
    import casadi.*
    
    %% initialize variables and jacobian/hessian functions
    
    all_m = 0;
    for i = 1:N
        lagrangian{i} = 0;
        cost{i} = 0;
        private_eq_constraints{i} = [];
        private_ineq_constraints{i} = [];
        control_vars{i} = [];
        all_m = all_m + m{i};
    end
    
    eps = MX.sym('eps', 1);
    cost{p} = 10*eps;
    x{1} = MX.sym(['x_' num2str(1)], n);
    
    state_vars = [x{1}];
    shared_dyn_constraints = [z0-x{1}];
    
    for t = 1:T
        x{t+1} = MX.sym(['x_' num2str(t+1)], n);
        state_vars = [state_vars; x{t+1}];
        
        xu = [x{t}];
        for i = 1:N
            u{t,i} = MX.sym(['u_' num2str(t) '_' num2str(i)], m{i});
            control_vars{i} = [control_vars{i}; u{t,i}];
            xu = [xu; u{t,i}];
        end
        for i = 1:N
            private_eq_constraints{i} = [private_eq_constraints{i}; h{t,i}(xu)];
            private_ineq_constraints{i} = [private_ineq_constraints{i}; g{t,i}(xu)];
        end
        
        for i = 1:N
            cost{i} = cost{i} + l{t,i}(xu);
        end
        
        shared_dyn_constraints = [shared_dyn_constraints; f(xu)-x{t+1}];  
    end
    
    for i = 1:N
        private_eq_constraints{i} = [private_eq_constraints{i}; h{T+1,i}(x{T+1})];
        private_ineq_constraints{i} = [private_ineq_constraints{i}; g{T+1,i}(x{T+1})];
        cost{i} = cost{i} + l{T+1,i}(x{T+1});
    end
    
    num_seq = size(shared_dyn_constraints,1);
    
    eq_constraints = [shared_dyn_constraints;
                      private_eq_constraints{p}];
    ineq_constraints = [private_ineq_constraints{p}; -eps];
    primal_vars = [state_vars; eps; control_vars{p}];
     
    for i = 1:N
        lagrangian{i} = cost{i};
        if i ~= p
            num_peq{i} = size(private_eq_constraints{i},1);
            num_pineq{i} = size(private_ineq_constraints{i},1);

            mu{i} = MX.sym(['mu_' num2str(i)],num_peq{i});
            lam{i} = MX.sym(['lam_' num2str(i)],num_seq);
            gam{i} = MX.sym(['gam_' num2str(i)],num_pineq{i});

            if num_peq{i} > 0
                lagrangian{i} = lagrangian{i} + mu{i}'*private_eq_constraints{i};
            end
            if num_pineq{i} > 0
                lagrangian{i} = lagrangian{i} + gam{i}'*private_ineq_constraints{i};
            end
            if num_seq > 0
                lagrangian{i} = lagrangian{i} + lam{i}'*shared_dyn_constraints;
            end
            
            eq_constraints = [eq_constraints; 
                              gradient(lagrangian{i},state_vars);
                              gradient(lagrangian{i},control_vars{i});
                              private_eq_constraints{i}];
            ineq_constraints = [ineq_constraints;
                                private_ineq_constraints{i};
                                -gam{i};
                                -gam{i}'*private_ineq_constraints{i}-eps];
                                
            primal_vars = [primal_vars; 
                           control_vars{i}; 
                           mu{i}; 
                           lam{i};
                           gam{i}];
        end
    end
    
    extract_state = Function('state', {primal_vars}, {state_vars});
    
    grad_eq_constraints = jacobian(eq_constraints, primal_vars);
    grad_ineq_constraints = jacobian(ineq_constraints, primal_vars);
    eval_eq = Function('eq', {primal_vars}, {eq_constraints});
    eval_ineq = Function('ineq', {primal_vars}, {ineq_constraints});
    eval_grad_eq = Function('grad_eq', {primal_vars}, {grad_eq_constraints});
    eval_grad_ineq = Function('grad_ineq', {primal_vars}, {grad_ineq_constraints});
    
    lam{p} = MX.sym('mu', size(eq_constraints,1));
    gam{p} = MX.sym('gam', size(ineq_constraints,1));
    
    all_vars = [primal_vars; lam{p}; gam{p}];
    
    lagrangian{p} = cost{p} + lam{p}'*eq_constraints + gam{p}'*ineq_constraints;
    grad_primal_lagrangian = gradient(lagrangian{p}, primal_vars);
    hess_lag = jacobian(grad_primal_lagrangian, primal_vars);
    eval_hess = Function('hess_lag', {primal_vars, lam{p}, gam{p}}, {hess_lag});
    
    grad_cost = gradient(cost{p},primal_vars);
    eval_cost = Function('cost', {primal_vars}, {cost{p}});
    eval_grad_cost = Function('gradcost', {primal_vars}, {grad_cost});
    
    num_hess = @(x,l) sparse(eval_hess(x,l.eqnonlin,l.ineqnonlin));
    num_obj = @(x) helper_eval(x, eval_cost, eval_grad_cost);
    num_con = @(x) con_helper_eval(x, eval_ineq, eval_eq, eval_grad_ineq, eval_grad_eq);
    
    options = optimoptions('fmincon','Algorithm','interior-point',...
    "SpecifyConstraintGradient",true,"SpecifyObjectiveGradient",true,...
    'HessianFcn',num_hess);

    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    if length(varargin) == 0
        x0 = zeros(size(primal_vars));
        ind = 0;
        x0(ind+1:ind+n) = z0;
        ind = ind+n;
        for t = 1:T
            x0(ind+1:ind+n) = f([x0(ind-n+1:ind);0.0*randn(all_m,1)]);
            ind = ind+n;
        end
    else
        x0 = zeros(size(primal_vars));
        x0(1:size(varargin{1},1)) = varargin{1};
    end
    x0(size(state_vars,1)+1) = 5;
    
    tic;
    [x, fval, eflag, output] = fmincon(num_obj,x0,A,b,Aeq,beq,lb,ub,num_con,options);
    duration = toc;
    
    num_x_opt = full(extract_state(x));
    ind = 0;
    for t = 1:T+1
        x_opt(1:n,t) = num_x_opt(ind+1:ind+n);
        ind = ind+n;
    end
    eps_opt = x(ind+1);
    disp('done');

end