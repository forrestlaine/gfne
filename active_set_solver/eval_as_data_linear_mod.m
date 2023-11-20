function [F,H,G,Q] = eval_as_data_linear_mod(evaluators, xval,uval,lamval,muval,gamval,T,N,F,H,Q)
    for t = 1:T
        uu = [];
        for i = 1:N
            uu = [uu; uval{t,i}];
        end
        xu = [xval{t};uu];
        for i = 1:N
            H{t,i}(:,1) = 0;
            G{t,i}(:,1) = 0;
        end
        
        F{t}(:,1) = 0;

        for i = 1:N
            
            grad = full(evaluators.eval_grad_cost{t,i}(xu));
            hess = evaluators.eval_hess_cost{t,i}(xu);
            hess = full(hess);
            
            Q{t,i} = [0 grad'; 
                      grad hess];
        end
    end
    for i = 1:N
        H{T+1,i}(:,1) = 0;
        G{T+1,i}(:,1) = 0;
        
        grad = full(evaluators.eval_grad_cost{T+1,i}(xval{T+1}));
        hess = evaluators.eval_hess_cost{T+1,i}(xval{T+1});
        hess = full(hess);
        Q{T+1,i} = [0 grad'; 
                    grad hess];
    end  
end