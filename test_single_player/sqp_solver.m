%% Iterative EC solver

function [residuals, xval,uval] = sqp_solver(f,... % dynamics function
                                        e,... % shared ineq constraint function
                                        h,... % time-indexed dict of player-indexed equality constraint functions
                                        g,... % time-indexed dict of player-indexed inequality constraint functions
                                        l,... % time-indexed dict of player-indexed cost functions
                                        n,... % dimension of shared state
                                        m,... % player-indexed dict of player input dimensions
                                        N,... % number of players
                                        T,... % time steps
                                        rr,... % flags for non-dynamic plalyers
                                        z0)   % initial state
    import casadi.*
    %% initialize variables and jacobian/hessian functions

    all_m = 0;
    for i = 1:N
        lagrangian{i} = 0;
        private_eq_constraints{i} = [];
        private_ineq_constraints{i} = [];
        control_vars{i} = [];
        all_m = all_m + m{i};
    end
    
    x{1} = MX.sym(['x_' num2str(1)], n);
    ic_negative_slacks = MX.sym('ic_negative_slacks',n);
    ic_positive_slacks = MX.sym('ic_positive_slacks',n);
    
    state_vars = [x{1}];
    shared_dyn_constraints = [z0-x{1}];
    shared_ineq_constraints = [];
    
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
            lagrangian{i} = lagrangian{i} + l{t,i}(xu);
        end
        
        shared_dyn_constraints = [shared_dyn_constraints; f(xu)-x{t+1}];
        shared_ineq_constraints = [shared_ineq_constraints; e(x{t+1})];

       
    end
    
    for i = 1:N
        private_eq_constraints{i} = [private_eq_constraints{i}; h{T+1,i}(x{T+1})];
        private_ineq_constraints{i} = [private_ineq_constraints{i}; g{T+1,i}(x{T+1})];
        lagrangian{i} = lagrangian{i} + l{T+1,i}(x{T+1});
    end
    
    
    num_seq = size(shared_dyn_constraints,1);
    num_sineq = size(shared_ineq_constraints,1);
    
    
    eta = MX.sym('eta]',num_sineq);
    
    for i = 1:N
        num_peq{i} = size(private_eq_constraints{i},1);
        num_pineq{i} = size(private_ineq_constraints{i},1);

        mu{i} = MX.sym(['mu_' num2str(i)],num_peq{i});
        gam{i} = MX.sym(['gam_' num2str(i)],num_pineq{i});


        if num_peq{i} > 0
            lagrangian{i} = lagrangian{i} - mu{i}'*private_eq_constraints{i};
        end
        if num_pineq{i} > 0
            lagrangian{i} = lagrangian{i} - gam{i}'*private_ineq_constraints{i};
        end
        if ~rr{i}

            lam{i} = MX.sym(['lam_' num2str(i)],num_seq);
            if num_seq > 0
                lagrangian{i} = lagrangian{i} - lam{i}'*shared_dyn_constraints;
            end
            if num_sineq > 0
                lagrangian{i} = lagrangian{i} - eta'*shared_ineq_constraints;
            end
        end
    end
    
    foncs = [shared_dyn_constraints];
    all_vars = [state_vars];
    lower_bound = [-ones(num_seq,1)];

    for i = 1:N
        if ~rr{i}
            foncs = [foncs; 
                    gradient(lagrangian{i},state_vars);
                    gradient(lagrangian{i},control_vars{i});
                    private_eq_constraints{i}];
            all_vars = [all_vars;
                        lam{i};
                        control_vars{i};
                        mu{i}];
            lower_bound = [lower_bound;
                           -ones(num_seq,1);
                           -ones(size(control_vars{i}));
                           -ones(num_peq{i},1)];
        else
            foncs = [foncs; 
                    gradient(lagrangian{i},control_vars{i});
                    private_eq_constraints{i}];
            all_vars = [all_vars;
                        control_vars{i};
                        mu{i}];
            lower_bound = [lower_bound;
                           -ones(size(control_vars{i}));
                           -ones(num_peq{i},1)];
        end
    end
    
    % number of variables not in complementarity conditions
    size_free_vars = size(all_vars,1); 
    
    foncs = [foncs; shared_ineq_constraints];
    all_vars = [all_vars; eta];
    lower_bound = [lower_bound; zeros(num_sineq,1)];
    all_ineq_mults = [eta];
    
    for i = 1:N
        foncs = [foncs; private_ineq_constraints{i}];
        all_vars = [all_vars; gam{i}];
        all_ineq_mults = [all_ineq_mults; gam{i}];
        lower_bound = [lower_bound; zeros(num_pineq{i},1)];
    end
    
    size_all_vars = size(all_vars,1);
    for i = 1:4
        extract_control_vars{i} = Function('extract_controls', {all_vars}, {control_vars{i}});
    end
             
    % Complementarity system is set up now. 
    % MCP([lower_bound, upper_bound], foncs)
    % MCP([a,b],F) <-> {     (Fi(z) = 0 and ai <= zi <= bi) 
    %                     or (Fi(z) > 0 and zi = ai)
    %                     or (Fi(z) < 0 and zi = bi)          }
    % lower_bound looks like [-inf ... -inf 0 ... 0]
    % upper_bound looks like [inf      ...      inf]
    
    H = jacobian(foncs, all_vars);
    M_free_vars = H(1:size_free_vars,1:size_free_vars);
    m_free_vars = foncs(1:size_free_vars);
    N_free_vars = H(1:size_free_vars,size_free_vars+1:end);
    G = H(size_free_vars+1:end,1:size_free_vars);
    
    QP_hess = -G/M_free_vars*N_free_vars;
    QP_grad = -G/M_free_vars*m_free_vars + foncs(size_free_vars+1:end);

    eval_F = Function('eval_foncs', {all_vars}, {foncs});
    eval_H = Function('eval_hess', {all_vars}, {H});
    
    eval_QP_hess = Function('eval_QP_hess', {all_vars}, {QP_hess});
    eval_QP_grad = Function('eval_QP_grad', {all_vars}, {QP_grad});
    
    
    constrained_vars = MX.sym('constrained_vars', size_all_vars-size_free_vars);
    free_vars = -M_free_vars\(N_free_vars*constrained_vars + m_free_vars);
    eval_free_vars = Function('eval_free_vars', {constrained_vars, all_vars}, {free_vars});
    
    z = 0*randn(size_all_vars,1);
    x0 = [z0];
    for t = 1:T
        x0 = [x0; f([x0(end-n+1:end);10*randn(all_m,1)])];
    end
    z(1:size(state_vars,1)) = x0;

    resid = inf;
    
    options =  optimoptions('quadprog');
%     options =  optimoptions('quadprog','Display','off');
%     options =  optimoptions('quadprog','OptimalityTolerance',1e-3);

    for t = 1:100
        if size(constrained_vars,1) > 0
%             zhat = z;
            z(size_free_vars+1:end) = 0;
%             HH = full(eval_H(z));
%             MM = HH(1:size_free_vars,1:size_free_vars);
%             hh = full(eval_F(z));
%             mm = hh(1:size_free_vars);
%             gg = hh(size_free_vars+1:end);
%             NN = HH(1:size_free_vars,size_free_vars+1:end);
%             GG = HH(size_free_vars+1:end,1:size_free_vars);
            
%             Q = -GG*lsqminnorm(MM,NN);
%             q = -GG*lsqminnorm(MM,mm) + gg;
            
            Q = full(eval_QP_hess(z));
            q = full(eval_QP_grad(z));
            
            num_g = size(q,1);
            
%             H = blkdiag(Q+Q', 0*eye(1*num_g));
            QQ = Q+Q';
            QQ = Q+Q' + (abs(min(real(eig(QQ))))+0.0001)*eye(size(QQ,1));
            H = sparse([QQ zeros(num_g); zeros(num_g) zeros(num_g)]);
            h = [q; 1000000*ones(num_g,1)];
%             min_eig = min(real(eig(H)));
            
        
            A = sparse([-Q -eye(num_g)]);
            b = [q];
            
            lb = [zeros(num_g,1); zeros(num_g,1)];
            ub = Inf*ones(size(lb));
            
            
            disp('Solving QP:');
            model.Q = sparse(H/2);
            model.obj = h;
            model.A = sparse(A);
            model.rhs = b;
            model.sense = '<';
            model.lb = lb;
            model.ub = ub;
            params.OutputFlag = 0;
            
            results = gurobi(model,params);
            
            
            
%             [delta_gam,fval,flag] = quadprog(Q+Q',q,-Q,q,[],[],zeros(size(q)),Inf*ones(size(q)),[],options);
%             [delta_gam_slacks,fval,flag,output,lambdas] = quadprog(H,h,A,b,[],[],lb,ub,[],options);
%             if flag ~= 1 || fval > 1e-1
%                 % flag  |  status
%                 %   1       Good
%                 %  -2       Infeasible
%                 %  -3       Unbounded
%                 %  -6       Non-convex
% %                 min_eig
%                 disp('QP NOT FEASIBLE!');
%                 fval
%             end
            delta_gam_slacks = results.x;
            delta_gam = delta_gam_slacks(1:num_g);
            slacks = delta_gam_slacks(num_g+1:end);
%             delta_free_vars = -lsqminnorm(MM,NN*delta_gam + mm);
            delta_free_vars = eval_free_vars(delta_gam,z);
            delta = [delta_free_vars; delta_gam];
        else
            delta = full(-eval_H(z)\eval_F(z));
        end
        
        if t==1
            alpha = 1;
        else
            alpha = 1;
        end
        
        linesearch_suceeded = resid > 1e3;
        disp('need to achieve');
        resid
        for s = 1:10 % max linesearch iters
            candidate = z + alpha*delta;
%             candidate(size_free_vars+1:end) = max(candidate(size_free_vars+1:end),0);
            candidate_foncs = full(eval_F(candidate));
            candidate_foncs(size_free_vars+1:end) = min(candidate_foncs(size_free_vars+1:end),0);
%             candidate_resid = max(norm(candidate_foncs(1:size_free_vars)), norm(candidate_foncs(size_free_vars+1:end)));
            candidate_resid = norm(candidate_foncs,1);
            if candidate_resid < resid
                linesearch_suceeded = true;
                break;
            else
                candidate_resid
                alpha = alpha * 0.5;
            end
        end
        
        if linesearch_suceeded
            resid = candidate_resid;
            z = candidate;
            actual_foncs = candidate_foncs;
        else
            disp('Linesearch did not converge!');
            break;
        end
        
        resid
        if resid < 1e-5
            break;
        end
    end
    
    for i = 1:N
        U_num{i} = extract_control_vars{i}(z);
        U{i} = zeros(m{i},T);
        indi{i} = 0;
    end
    
    X = zeros(n,T+1);
    
    
    ind = 0;
    
    for t=1:T
        X(:,t) = full(z(ind+1:ind+n));
        ind = ind+n;
        for i = 1:N
            U{i}(:,t) = full(U_num{i}(indi{i}+1:indi{i}+m{i}));
            indi{i} = indi{i}+m{i};
        end
    end
    X(:,T+1) = full(z(ind+1:ind+n));
   
    xval = X;
    residuals = resid;
    uval = U;
    
end