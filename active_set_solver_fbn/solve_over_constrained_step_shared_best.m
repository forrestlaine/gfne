function [G,GI,P,upol] = solve_over_constrained_step_shared_best(Ft,...
                                                  H,...
                                                  HI,...
                                                  GG,...
                                                  GI,...
                                                  Qt,... 
                                                  P,...
                                                  N,... 
                                                  m,...
                                                  full_scope)   
      
    n = size(Ft,1);                                         
    all_m = 0;
    for i = 1:N
        m_start{i} = all_m;
        all_m = all_m + m{i};
    end

    dyn = [1 zeros(1,n+all_m);
           Ft];
    
       
    Hall = [H; GG*dyn];
    HIall = [HI, GI];
       
    Mall = [];
    Rall = [];
    RIall = {};
    RIiter = 1;
    tot = 0;
    
    for i = 1:N
        Hti = zeros(0,n+1+all_m);
        Gi = zeros(0,n+1);
        for j = 1:size(H,1)
            if ismember(i,HI{j})
                Hti = [Hti;H(j,:)];
            end
        end
        for j = 1:size(GG,1)
            if ismember(i,GI{j})
                Gi = [Gi; GG(j,:)];
            end
        end
        Ht{i} = Hti;
        G{i} = Gi;
    end
    
    for i = 1:N
        Qto{i} = Qt{i};
        Qt{i} = Qt{i} + dyn'*P{i}*dyn;
        nh = size(Ht{i},1);
        Ht{i} = [Ht{i}; G{i}*dyn];
        Hut{i} = Ht{i}(:,1+n+m_start{i}+1:1+n+m_start{i}+m{i});
        [U,S,V] = svd(Hut{i});
        rk{i} = rank(Hut{i},1e-6);
        z{i} = m{i}-rk{i};
        n_start{i} = tot;
        tot = tot+z{i};
        Ur = U(:,1:rk{i});
        Un = U(:,rk{i}+1:end);
        Vr{i} = V(:,1:rk{i});
        Vn{i} = V(:,rk{i}+1:end);
        Sr = S(1:rk{i},1:rk{i});
        Mt{i} = Ur'*Ht{i};
        Mall = [Mall; Mt{i}];
        Rt{i} = Un'*Ht{i};
        Rall = [Rall; Rt{i}];
        for j = 1:size(Rt{i},1)
            RIall{RIiter} = i;
            RIiter = RIiter + 1;
        end
    end

    Vra = blkdiag(Vr{:});
    Vna = blkdiag(Vn{:});
    Va = [Vra Vna];

    mod = [eye(1+n) zeros(1+n,all_m);
           zeros(all_m, 1+n) Va];
    % [1;x;u] = mod*[1;x;z1;z2]
    Mall = Mall*mod;

    if rank(Mall(:,1+n+1:1+n+size(Vra,2)), 1e-5) < size(Mall(:,1+n+1:1+n+size(Vra,2)),1)
        disp('Rank degenerate in over_constrained (Incompatable constraints)');
    end
    z1pol = -Mall(:,1+n+1:1+n+size(Vra,2))\ [Mall(:,1:1+n) Mall(:,2+n+size(Vra,2):end)];
    % z1 = z1pol * [1;x;z2]

    z2mod = mod*[eye(1+n) zeros(1+n, size(z1pol,2)-1-n);
                 z1pol;
                 zeros(size(z1pol,2)-1-n,1+n), eye(size(z1pol,2)-1-n)];
    % [1;x;u] = z2mod*[1;x;z2]

    Nall = [];
    for i = 1:N
        Qtz{i} = z2mod'*Qt{i}*z2mod;
        Nall = [Nall; Qtz{i}(1+n+n_start{i}+1:1+n+n_start{i}+z{i},:)];
    end
    if rank(Nall(:,1+n+1:end),1e-5) < size(Nall(:,1+n+1:end),1)
        disp('Rank degenerate in over_constrained (Incompatable objectives)');
    end
    
    z2pol = -Nall(:,1+n+1:end)\Nall(:,1:1+n);
    % z2 = z2pol * [1;x]

    upol = Va*[z1pol;
               zeros(size(z2pol,1),1+n) eye(size(z2pol,1))] * [eye(n+1);z2pol];

    % u = upol*[1;x]
    upd = [eye(1+n);upol];
    
%     allH = [H;GG*dyn]*upd;
    allH = Rall*upd;
    allHI = RIall;
%     for j = 1:size(allH,1)
%         if norm(allH(j,:)) > 1e-6
%             allH(j,:) = allH(j,:) / norm(allH(j,:));
%         end
%     end
           
%     allHI = [HI, GI];
    
    for i = 1:N
        if full_scope{i}
            P{i} = upd'*Qt{i}*upd;
            P{i} = 0.5*(P{i} + P{i}');
        else
            P{i} = Qto{i}(1:n+1,1:n+1);
        end
        P{i}(1,1) = 0;
    end
    
    G = zeros(0,n+1);
    GI = {};
    
%     j = 1;
%     for i = 1:N
%         G = [G; Rt{i}*upd];
%         for l = 1:size(Rt{i},1)
%             GI{j} = i;
%             j = j+1;
%         end
%     end
            
    i = 1;
    for j = 1:size(allH,1)
        if rank(allH(j,:),1e-6) > 0
            G = [G; allH(j,:)];
            GI{i} = allHI{j};
            i = i+1;
        end
    end
end