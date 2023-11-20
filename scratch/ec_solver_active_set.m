%% Iterative EC solver

function [residuals, xval,full_sol] = ec_solver_active_set(f,... % dynamics function
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
  
    all_m = 0;
    for i = 1:N
        start_m{i} = all_m;
        all_m = all_m + m{i};
    end
    
    evaluators = generate_evaluators(f,h,g,l,n,m,N,T);
    [xval,...
     uval,...
     lamval,...
     muval,...
     slackval,...
     gamval,...
     devval,...
     devopt,...
     psival,...
     tauval] = initialize_sol(params,z0,T,N,n,m,all_m,evaluators,initialization);
 
    [F,H,G,Q] = eval_as_data(evaluators,xval,uval,lamval,muval,gamval,T,N);
    H_active = objcopy(H);
    [K,k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
    [XX,UU,LL,MM_active,PP] = solve_ec_lq_game_d(F,H,Q,N,T,K,k,zeros(n,1));
    for t = 1:T
        for i = 1:N
            working_set{t,i} = false(size(gamval{t,i})); 
        end
    end
    
    while true
        while true
            for t = 1:T+1
                if t < T+1
                    vec = [1;XX{t};UU{t,i}];
                else
                    vec = [1;XX{t}];
                end
                for i = 1:N
                    for j = 1:size(gamval{t,i},1)
                        if G(j,:)*vec < 0
                            working_set{t,i}(j) = true;
                            H_active{t,i} = [H_active{t,i}; G(j,:)];
                        end
                    end
                end
            end   
            [K,k] = solve_ec_lq_game_r(F,H_active,Q,N,m,T);
            [dXX,dUU,LL,MM_active,PP] = solve_ec_lq_game_d(F,H_active,Q,N,T,K,k,zeros(n,1));
            step = 0;
            for t = 1:T+1
                step = step + norm(dXX{t});
            end
            if step < 1e-3
                break;
            end
        end
    end
    
   
   for iter = 1:100
       
        [active_F, active_H, active_Q,xval,uval
        for t = 1:T
            for i = 1:N
                kd{t,i} = K{t,i}*xval{t}+k{t,i}+devval{t,i}-uval{t,i};
                devdelt{t,i} = devval{t,i}-devopt{t,i};
            end
        end
        [dX,dU,dL,dM,dP,dD] = solve_ec_lq_game_regdev(F,H,Q,N,T,K,kd,zeros(size(z0)),0.0,devval);

        solve_time = toc;
        
        tic;
        [alpha_dual, alpha_primal, dG, dS] = extract_slacks_and_gams(params,dX,dU,xval,uval,gamval,slackval,tauval,evaluators,T,N);
        extract_slack_time = toc;

        deltas.dX = dX;
        deltas.dU = dU;
        deltas.dL = dL;
        deltas.dM = dM;
        deltas.dP = dP;
        deltas.dD = dD;
        deltas.dG = dG;
        deltas.dS = dS;
        %% linesearch (update this so we're not packing and upacking data)
       
        if iter <= params.warm_up_iters
            alpha = params.warm_up_alpha;
        else
            alpha = params.alpha_init;
        end
        start_s = 1;

        if params.debug_print >= 2
            alpha_primal
            alpha_dual
        end

        linesearch_succeeded = false;
        if params.open_loop || ~ol_converged
            num_linesearch_iters = params.max_linesearch_iters;
        else
            num_linesearch_iters = params.max_cl_linesearch_iters;
        end
        
        tic;
        for s = start_s:num_linesearch_iters % max linesearch iterations
            % pack candidate values ('c_') and jacobians, hessians, etc.
            [c_xval,...
            c_uval,...
            c_devval,...
            c_lamval,...
            c_muval,...
            c_psival,...
            c_slackval,...
            c_gamval] = extract_search_point(xval,...
                                            uval,...
                                            devval,...
                                            lamval,...
                                            muval,...
                                            psival,...
                                            slackval,...
                                            gamval,...
                                            deltas,...
                                            alpha,...
                                            alpha_primal,...
                                            alpha_dual,...
                                            params,...
                                            restoration,...
                                            evaluators,...
                                            T, N);
            new_zvec = pack_zvec(c_xval,c_uval,c_muval,c_lamval,c_gamval,c_slackval,T,N);
            [c_current_residual,...
             c_private_residual,...
             c_private_resids,...
             c_constraint_residual] = evaluate_step(new_zvec,...
                                                    evaluators,...
                                                    N,T,n,m,all_m,start_m,...
                                                    K,k,...
                                                    tauval,...
                                                    c_xval,...
                                                    c_uval,...
                                                    c_devval,...
                                                    devopt,...
                                                    c_psival,...
                                                    c_gamval,...
                                                    c_slackval,...
                                                    reg);

           
            if (c_current_residual < current_residual)
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

  
        current_residual = c_current_residual;
        private_residual = c_private_residual;
        private_fine_resids = c_private_resids;
        all_constraint_residual = c_constraint_residual;
        
        xval = c_xval;
        uval = c_uval;
        devval = c_devval;
        lamval = c_lamval;
        muval = c_muval;
        psival = c_psival;
        gamval = c_gamval;
        slackval = c_slackval;
        zvec = new_zvec;
        
        if params.debug_print >= 1
            current_residual
            tauval
        end
        
        [min_comp, comp, dim] = eval_complementarity_stats(N,T,gamval,slackval);
        total_dim = 0;
        for i = 1:N
            total_dim = total_dim + dim{i};
        end

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
        all_taus_down = all_taus_down || total_dim == 0;
        
        constraints_satisfied = all_constraint_residual < params.constraint_tolerance;
        
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

                [F,H,Q] = eval_data(evaluators, xval,uval,lamval,muval,slackval,gamval,tauval,T,N);
                [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
%                 for t = 1:T
%                     for i = 1:N
%                         if slackval{t,i} < 1e-3
%                             for tt = t:-1:max(1,t-3)
%                                 for j = 1:N
%                                     K{tt,j} = zeros(size(K{tt,j}));
%                                     k{tt,j} = zeros(size(k{tt,j}));
%                                 end
%                             end
%                         end
%                     end
%                 end
                fresh_policies = true
                fresh_policy_count = fresh_policy_count + 1;
                reg = params.dev_penalty * params.dev_decay^(fresh_policy_count-1);

%                 for i = 1:N
%                     tauval{i} = params.tauval_init{i};
%                 end
                for t = 1:T
                    for i = 1:N
                        devopt{t,i} = uval{t,i}-(K{t,i}*xval{t}+k{t,i});
                        devval{t,i} = devopt{t,i};
                        psival{t,i}(start_m{i}+1:start_m{i}+m{i}) = -(devval{t,i}-devopt{t,i})*params.dev_penalty;
                    end
                end
                
             [current_residual,...
              private_residual,...
              private_resids,...
              constraint_residual] = evaluate_step(zvec,...
                                                   evaluators,...
                                                    N,T,n,m,all_m,start_m,...
                                                    K,k,...
                                                    tauval,...
                                                    xval,...
                                                    uval,...
                                                    devval,...
                                                    devopt,...
                                                    psival,...
                                                    gamval,...
                                                    slackval,...
                                                    reg);
                continue;
            elseif all_taus_down
                break;
            else
                tauval = update_tauval(params,restoration,min_comp,comp,dim,tauval,fresh_policy_count,current_residual,N);
            end
        else
            fresh_policies = false;
            tauval = update_tauval(params,restoration,min_comp,comp,dim,tauval,fresh_policy_count,current_residual,N);
            
            [F,H,Q] = eval_data(evaluators, xval,uval,lamval,muval,slackval,gamval,tauval,T,N);
            if ~params.open_loop && (~params.converge_on_every_policy || ol_converged)
                [K, k] = solve_ec_lq_game_r(F,H,Q,N,m,T);
                if params.converge_on_every_policy
                    ol_converged = false;
                end
                fresh_policies = true;
                fresh_policy_count = fresh_policy_count + 1;
            end
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