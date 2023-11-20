function [F,H,G,Q] = eval_as_data_linear(evaluators, xval,uval,lamval,muval,gamval,T,N,F,H,Q)
    for t = 1:T
        uu = [];
        for i = 1:N
            H{t,i}(:,1) = 0;
            G{t,i} = [full(evaluators.eval_ineq_constraint{t,i}([xval{t};uval{t,i}])),...
                      full(evaluators.eval_ineq_jac_constraint{t,i}([xval{t};uval{t,i}]))];
            uu = [uu; uval{t,i}];
        end
        xu = [xval{t};uu];
        F{t}(:,1) = 0;

        for i = 1:N
            
            grad = full(evaluators.eval_grad_cost{t,i}(xu));
            hess = evaluators.eval_hess_cost{t,i}(xu)+...
                   evaluators.eval_hess_ineq_constraint_mult{t,i}(xu,gamval{t,i});
            hess = full(hess);
            
            Q{t,i} = [0 grad'; 
                      grad hess];
        end
    end
    for i = 1:N
        H{T+1,i}(:,1) = 0;
        G{T+1,i} = [full(evaluators.eval_ineq_constraint{T+1,i}(xval{T+1})), full(evaluators.eval_jac_ineq_constraint{T+1,i}(xval{T+1}))];
        
        grad = full(evaluators.eval_grad_cost{T+1,i}(xval{T+1}));
        hess = evaluators.eval_hess_cost{T+1,i}(xval{T+1})+...
            evaluators.eval_hess_ineq_constraint_mult{T+1,i}(xval{T+1},gamval{T+1,i});
        hess = full(hess);
        Q{T+1,i} = [0 grad'; 
                    grad hess];
    end  
end