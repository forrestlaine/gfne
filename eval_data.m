function [F,H,Q] = eval_data(evaluators, xval,uval,lamval,muval,slackval,gamval,tauval,T,N)
    for t = 1:T
        uu = [];
        for i = 1:N
            H{t,i} = [full(evaluators.eval_constraint{t,i}([xval{t};uval{t,i}])),...
                      full(evaluators.eval_jac_constraint{t,i}([xval{t};uval{t,i}]))];
            uu = [uu; uval{t,i}];
        end
        xu = [xval{t};uu];
        F{t} = [full(evaluators.eval_pred{t}(xu))-xval{t+1}, full(evaluators.eval_jac_pred{t}(xu))];

        for i = 1:N
            full_ineq_jac{t,i} = evaluators.eval_full_jac_ineq_constraint{t,i}(xu);
            Sig = sparse(diag(gamval{t,i}./slackval{t,i}));
            ineq_const_val{t,i} = evaluators.eval_ineq_constraint{t,i}([xval{t};uval{t,i}]);
            grad = full(evaluators.eval_grad_cost{t,i}(xu) +...
                        full_ineq_jac{t,i}'*(Sig*ineq_const_val{t,i} -...
                                        gamval{t,i} -...
                                        tauval{i} * slackval{t,i}.^(-1)));
            hess = evaluators.eval_hess_cost{t,i}(xu)+...
                   evaluators.eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                   evaluators.eval_hess_dyn_mult{t,i}(xu,lamval{t,i})-...
                   evaluators.eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t,i});
            aug = full_ineq_jac{t,i}'*Sig*full_ineq_jac{t,i};
            hess = full(hess+aug);
            
            Q{t,i} = [0 grad'; 
                      grad hess];
        end
    end
    for i = 1:N
        H{T+1,i} = [full(evaluators.eval_constraint{T+1,i}(xval{T+1})), full(evaluators.eval_jac_constraint{T+1,i}(xval{T+1}))];
        
        full_ineq_jac{T+1,i} = evaluators.eval_jac_ineq_constraint{T+1,i}(xval{T+1});
        Sig = sparse(diag(gamval{T+1,i}./slackval{T+1,i}));
        ineq_const_val{T+1,i} = evaluators.eval_ineq_constraint{T+1,i}(xval{T+1});
        grad = full(evaluators.eval_grad_cost{T+1,i}(xval{T+1}) +...
                        full_ineq_jac{T+1,i}'*(Sig*ineq_const_val{T+1,i} -...
                                        gamval{T+1,i} -...
                                        tauval{i} * slackval{T+1,i}.^(-1)));
            
            
        hess = evaluators.eval_hess_cost{T+1,i}(xval{T+1})+...
            evaluators.eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i})-...
            evaluators.eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1,i});
        aug = full_ineq_jac{T+1,i}'*Sig*full_ineq_jac{T+1,i};
        hess = full(hess+aug);
        Q{T+1,i} = [0 grad'; 
                    grad hess];
    end  
end