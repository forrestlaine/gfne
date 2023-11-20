function [F,H,G,GR,Q] = eval_as_data_mod(evaluators, xval,uval,lamval,muval,gamval,T,N)
    for t = 1:T
        uu = [];
        for i = 1:N
            uu = [uu; uval{t,i}];
        end
        xu = [xval{t};uu];
        for i = 1:N
            H{t,i} = [full(evaluators.eval_constraint{t,i}(xu)),...
                      full(evaluators.eval_jac_constraint{t,i}(xu))];
            G{t,i} = [full(evaluators.eval_ineq_constraint{t,i}(xu)),...
                      full(evaluators.eval_ineq_jac_constraint{t,i}(xu))];
            GR{t,i} = [full(evaluators.eval_region{t,i}(xu)),...
                       full(evaluators.eval_jac_region{t,i}(xu))];
        end
        
        F{t} = [full(evaluators.eval_pred{t}(xu))-xval{t+1}, full(evaluators.eval_jac_pred{t}(xu))];

        for i = 1:N
            cost = full(evaluators.eval_cost{t,i}(xu));
            grad = full(evaluators.eval_grad_cost{t,i}(xu));
            hess = evaluators.eval_hess_cost{t,i}(xu)+...
                   evaluators.eval_hess_constraint_mult{t,i}(xu,muval{t,i})+...
                   evaluators.eval_hess_dyn_mult{t,i}(xu,lamval{t,i})+...
                   evaluators.eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t,i});
            hess = full(hess);
            
            Q{t,i} = [cost grad'; 
                      grad hess];
        end
    end
    for i = 1:N
        H{T+1,i} = [full(evaluators.eval_constraint{T+1,i}(xval{T+1})), full(evaluators.eval_jac_constraint{T+1,i}(xval{T+1}))];
        G{T+1,i} = [full(evaluators.eval_ineq_constraint{T+1,i}(xval{T+1})), full(evaluators.eval_jac_ineq_constraint{T+1,i}(xval{T+1}))];
        GR{T+1,i} = [full(evaluators.eval_region{T+1,i}(xval{T+1})), full(evaluators.eval_jac_region{T+1,i}(xval{T+1}))];
        
        cost = full(evaluators.eval_cost{T+1,i}(xval{T+1}));
        grad = full(evaluators.eval_grad_cost{T+1,i}(xval{T+1}));
        hess = evaluators.eval_hess_cost{T+1,i}(xval{T+1})+...
            evaluators.eval_hess_constraint_mult{T+1,i}(xval{T+1},muval{T+1,i})+...
            evaluators.eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1,i});
        hess = full(hess);
        Q{T+1,i} = [cost grad'; 
                    grad hess];
    end  
end