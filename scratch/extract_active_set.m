function [H_active,active_constraints] = extract_active_set(dX,dU,dG,T,N,H,G,active_constraints)
    % Gotta figure out this logic
    for t = 1:T
        for i = 1:N
            grad_active = H{t,i}(2:end,1);
            hess_active = H{t,i}(2:end,2:end);
            ineq_constraint = G{t,i}(2:end,:)*[1;dX{t};dU{t,i}];
            for l = 1:size(ineq_constraint,1)
                if active_constraint{t,i}(l) && dG{t,i}(l) <= 0
                    active_constraint{t,i}(l) = false;
                    
                end
            end
        end
    end

end