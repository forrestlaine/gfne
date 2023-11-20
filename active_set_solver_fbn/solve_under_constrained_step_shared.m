function [P,K,L,Lampol] = solve_under_constrained_step_shared(Ft,...
                                                  Ht,...
                                                  HIt,...
                                                  G,...
                                                  GI,...
                                                  Qt,... 
                                                  P,...
                                                  N,... 
                                                  m,...
                                                  all_m,...
                                                  full_scope) 
                                              
            n = size(Ft,1);                                                   
            dyn = [1 zeros(1,n+all_m);
                   Ft];
            sm = 0;
            for i = 1:N
                start_m{i} = sm;
                sm = sm+m{i};
            end

            all_lg = size(G,1);
            all_lh = size(Ht,1);
            for i = 1:N
                Gi{i} = G;
                Hti{i} = Ht;
                for j = 1:size(G,1)
                    if ~ismember(i,GI{j})
                        Gi{i}(j,:) = zeros(1,n+1);
                    end
                end
                for j = 1:size(Ht,1)
                    if ~ismember(i,HIt{j})
                        Hti{i}(j,:) = zeros(1,n+1+all_m);
                    end
                end
            end
            
            
            M = zeros(all_m + all_lh + (N+1)*n + all_lg);
            Mx = zeros(size(M,1),n);
            Mc = zeros(size(M,1),1);

            ind = 0;
            ind2 = 0;
            ind3 = 0;
            for i = 1:N
                M(ind+1:ind+m{i},1:all_m) = Qt{i}(ind+n+2:ind+n+m{i}+1,2+n:end);
                M(ind+1:ind+m{i},all_m+1:all_m+all_lh) = Hti{i}(:,2+n+start_m{i}:1+n+start_m{i}+m{i})';
                M(ind+1:ind+m{i},all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i) = Ft(:,2+n+ind:1+n+ind+m{i})';
                Mx(ind+1:ind+m{i},:) = Qt{i}(ind+n+2:ind+n+m{i}+1,2:1+n);
                Mc(ind+1:ind+m{i}) = Qt{i}(ind+n+2:ind+n+m{i}+1,1);
                ind = ind+m{i};
            end
            ind2 = 0;
            M(ind+1:ind+all_lh,1:all_m) = Ht(:,2+n:end);
            Mx(ind+1:ind+all_lh,:) = Ht(:,2:1+n);
            Mc(ind+1:ind+all_lh) = Ht(:,1);
            ind = ind+all_lh;
            
            
            M(ind+1:ind+n,1:all_m) = Ft(:,2+n:end);
            Mx(ind+1:ind+n,:) = Ft(:,2:1+n);
            Mc(ind+1:ind+n) = Ft(:,1);
            M(ind+1:ind+n,all_m+all_lh+n*N+1:all_m+all_lh+n*N+n) = -eye(n);
            ind = ind+n;
            ind2 = all_m+all_lh+n*N;
            for i = 1:N
                M(ind+1:ind+n,ind2+1:ind2+n) = P{i}(2:n+1,2:n+1);
                M(ind+1:ind+n,ind2-(N+1-i)*n+1:ind2-(N-i)*n) = -eye(n);
                M(ind+1:ind+n,ind2+n+1:ind2+n+all_lg) = Gi{i}(:,2:end)';
                Mc(ind+1:ind+n) = P{i}(2:n+1,1);
                ind = ind+n;
            end

            M(ind+1:ind+all_lg,ind2+1:ind2+n) = G(:,2:end);
            Mc(ind+1:ind+all_lg) = G(:,1);
            ind = ind+all_lg;
            
%             t = independent_systems(tt);
            M = sparse(M);
            Mrhs = sparse([Mc Mx]);
            if rank(full(M)) < size(M,1)
                disp('Rank degenerate in under_constrained!');
            end
            policies = -M\Mrhs;
            K = policies(1:all_m,:);
            ind1 = 0;
            ind2 = 0;
            L = -[policies(all_m+1:all_m+all_lh,:);
                  policies(all_m+all_lh+(N+1)*n+1:all_m+all_lh+(N+1)*n+all_lg,:)];
            for i = 1:N
                Lampol{i} = policies(all_m+all_lh+1+(i-1)*n:all_m+all_lh+n*i,:);
            end
            pol = [eye(n+1); K];
            
            for i = 1:N
                if full_scope{i}
                    P{i} = pol'*(Qt{i} + dyn'*P{i}*dyn)*pol;
                    P{i} = 0.5*(P{i} + P{i}');
                else
                    P{i} = Qt{i}(1:n+1,1:n+1);
                end
                P{i}(1,1) = 0;
            end
           
                
end