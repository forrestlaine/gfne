function [c_current_residual,...
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
                                                reg)
   
    dynamic_conds = evaluators.eval_dynamics(new_zvec);
    all_constraints = [dynamic_conds];
    new_necessary_conditions = [dynamic_conds];
    all_conds = [dynamic_conds];
    state_opt_conds = [];
    for i = 1:N
        % necessary conditions for state and other_control
        % optimalities need to be corrected to account for
        % anticipated-intent constraints/multipliers.
        constraint_opt = evaluators.eval_constraint_violation{i}(new_zvec);
        ineq_constraint_opt = evaluators.eval_ineq_constraint_violation{i}(new_zvec);
        raw_ineq_constraint_opt = evaluators.eval_ineq_constraint_violation{i}(new_zvec);
        state_opt = evaluators.eval_state_optimality{i}(new_zvec);
        control_opt = evaluators.eval_control_optimality{i}(new_zvec);
        all_control_opt = evaluators.eval_full_control_optimality{i}(new_zvec);
        other_control_opt = evaluators.eval_other_control_optimality{i}(new_zvec);
        complementarity_opt = [];
        policy_opt = [];
        deviation_opt = [];

        ind = 0;
        ind2 = 0;

        for t = 1:T
            KK = [];
            kk = [];
            for j = 1:N
                KK = [KK; K{t,j}];
                kk = [kk; k{t,j}];
            end
            deviation_opt = [deviation_opt; 
                             reg*(c_devval{t,i}-devopt{t,i}) + c_psival{t,i}(start_m{i}+1:start_m{i}+m{i})];
%                    
            if t > 1
                state_opt(ind+1:ind+n) = state_opt(ind+1:ind+n) + KK'*c_psival{t,i};
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
                           min(full(raw_ineq_constraint_opt),0);
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
        c_private_residual{i} = norm(full(private_conds{i}));
    end
    c_current_residual = norm(full(new_necessary_conditions));
    c_constraint_residual = norm(full(all_constraints));
                 
         
end