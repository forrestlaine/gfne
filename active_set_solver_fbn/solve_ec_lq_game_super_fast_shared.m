%% Test EC LQ
% This version computes the explicit riccati backup, without considering
% the lagrange multipliers. First substitutes in dynamics, then each player
% minimizes norm of their own constraints (may be suboptimal with respect 
% to their objective if there are multiple ways of doing this).
function [dX,dU,dL,dM,dG,dP,flag] = solve_ec_lq_game_super_fast_shared(F,...
                                                                  Heq,...
                                                                  HIeq,...
                                                                  Hineq,...
                                                                  HIineq,...
                                                                  Q,...
                                                                  N,...
                                                                  T,...
                                                                  x0,...
                                                                  m,...
                                                                  full_scope,...
                                                                  open_loop)  
        
                                              
                                                  
    Kk = cell(T,1);
    Ll = cell(T+1,1);

    n = size(F{1},1);
    all_m = 0;
    for i = 1:N
        start_m{i} =  all_m;
        all_m = all_m + m{i};
    end
   
    for t = 1:T+1
        H{t} = [Heq{t};Hineq{t}];
        HI{t} = [HIeq{t}, HIineq{t}];
    end
    
    if open_loop
        
        for t = 1:T
            for i = 1:N
                KK{t,i} = zeros(m{i},n);
            end
        end

        [XX,UU,LAM,MM,PP] = solve_ec_lq_game_dpol_shared(F,H,HI,Q,N,T,KK,m,full_scope);
        flag = 1;
        dX{1} = x0;
        for t = 1:T
            for i = 1:N
                dU{t,i} = UU{t,i}*[1;x0];
                dL{t,i} = LAM{t,i}*[1;x0];
                if t < T
                    dP{t,i} = PP{t,i}*[1;x0];
                end
            end
            mugam = MM{t}*[1;x0];
            dM{t} = mugam(1:size(Heq{t},1));
            dG{t} = mugam(size(Heq{t},1)+1:end);
            dX{t+1} = XX{t+1}*[1;x0];
        end
        mugam = MM{T+1}*[1;x0];
        dM{T+1} = mugam(1:size(Heq{T+1},1));
        dG{T+1} = mugam(size(Heq{T+1},1)+1:end);
        
    else
    
        G{T} = H{T+1};
        GI{T} = HI{T+1};
        P{T} = {};
        for i = 1:N
            P{T}{i} = Q{T+1,i}; % cost-to-go
        end

        over_constrained = false;
        for t = T:-1:1
            all_feas{t} = true;
            % This check is over-restrictive but okay for now
            for i = 1:N
                B = F{t}(:,1+n+start_m{i}+1:1+n+start_m{i}+m{i});
                Hu = zeros(0,m{i});
                Gi = zeros(0,n+1);
                for j = 1:size(H{t},1)
                    if ismember(i,HI{t}{j})
                        Hu = [Hu; H{t}(j,n+2+start_m{i}:n+1+start_m{i}+m{i})];
                    end
                end
                if size(Hu,1) >= m{i}
                    ct = 0;
                    for jj = 1:numel(HI{t+1})
                        if ismember(i,HI{t+1}{jj})
                            ct = ct+1;
                        end
                    end
    %                 if ct > m{i}
    %                     disp('lets check this out');
    %                 end
                end
                for j = 1:size(G{t},1)
                    if ismember(i,GI{t}{j})
                        Gi = [Gi; G{t}(j,:)];
                    end
                end
                HHu = [Gi(:,2:end)*B; Hu];
                rk = rank(HHu);
                all_feas{t} = all_feas{t} && (rk == size(HHu,1));
            end

            if all_feas{t}
                if over_constrained
                    Ft = F{t};
                    dyn = [1 zeros(1,n+all_m);
                           Ft];
                    for i = 1:N
                        Qt{i} = Q{t,i};
                    end
                    Ht = H{t};
                    HIt = HI{t};

                    [~,~,Po,KKt] = solve_over_constrained_step_shared(Ft,Ht,HIt,G{t},GI{t},Qt,P{t},N,m,full_scope);
                    Kk{t} = KKt;

                    S = 1 + first_over_constrained - t;
                    for s = 1:S
                        FF{s} = F{t+s-1};
                        HH{s} = H{t+s-1};
                        HHI{s} = HI{t+s-1};
                        for i = 1:N
                            QQ{s,i} = Q{t+s-1,i};
                            KK{s,i} = Kk{t+s-1}(start_m{i}+1:start_m{i}+m{i},2:end);
                        end
                    end
                    HH{S+1} = Goc;
                    HHI{S+1} = GIoc;
                    for i = 1:N
                        QQ{S+1,i} = Poc{i};
                    end
                    [~,~,LAM,MM,~] = solve_ec_lq_game_dpol_shared(FF,HH,HHI,QQ,N,S,KK,m,full_scope);

                    if t+S == T+1
                        Ll{t+S} = MM{S+1};
                        Lid{t+S} = t;
                    end

                    for i = 1:N
                       LAMLAM{t+S-1,i} = LAM{S,i};
                    end

                    for s = S:-1:1
                        Ll{t+s-1} = MM{s};
                        Lid{t+s-1} = t;
                        for i = 1:N
                            if s>1
                                LAMLAM{t+s-2,i} = LAM{s-1,i};
                            end
                        end

                    end
                    if t > 1
                        P{t-1} = Po;
                        G{t-1} = zeros(0,n+1);
                        GI{t-1} = {};
                    end
                    % need to go and get multipliers
                else


                    Ft = F{t};
                    dyn = [1 zeros(1,n+all_m);
                           Ft];
                    Ht = H{t};
                    HIt = HI{t};
                    for i = 1:N
                        Qt{i} = Q{t,i};
                    end
                    if size(Ht,1) ~= numel(HIt)
                        disp('bad!');
                    end
                    [Po,KKt,LLt,Lampolt] = solve_under_constrained_step_shared(Ft,Ht,HIt,G{t},GI{t},Qt,P{t},N,m,all_m,full_scope);
                    Kk{t} = KKt;
                    if t > 1
                        G{t-1} = zeros(0,n+1);
                        GI{t-1} = {};
                        P{t-1} = Po;
                    end
                    Ll{t} = LLt(1:size(Ht,1),:);
                    if t==T
                        Ll{t+1} = LLt(size(Ht,1)+1:end,:);
                        Lid{t+1} = t;
                    end
                    for i = 1:N
                        LAMLAM{t,i} = Lampolt{i};
                    end
                    Lid{t} = t;

                end
                over_constrained = false;
            else
                if ~over_constrained
                    first_over_constrained = t;
                    Goc = G{t};
                    GIoc = GI{t};
                    Poc = P{t};
                end
                over_constrained = true;
                Ft = F{t};
                dyn = [1 zeros(1,n+all_m);
                       Ft];
                Ht = H{t};
                HIt = HI{t};
                for i = 1:N
                    Qt{i} = Q{t,i};
                end
                [Go,GIo,Po,KKt] = solve_over_constrained_step_shared(Ft,Ht,HIt,G{t},GI{t},Qt,P{t},N,m,full_scope);
                Kk{t} = KKt;
                if t > 1
                    G{t-1} = Go;
                    GI{t-1} = GIo;
                    P{t-1} = Po;
                end
            end
        end
        if ~exist('Lid','var')
            disp('Looks like infeasible problem');
            flag = 0;
            dX = 0 ;
            dU = 0;
            dL = 0;
            dM = 0;
            dG = 0;
            dP = 0;
        else
            dX{1} = x0;
            for t = 1:T
                u = Kk{t}*[1;dX{t}];
                mugam = Ll{t}*[1;dX{Lid{t}}];
                dM{t} = mugam(1:size(Heq{t},1),:);
                dG{t} = mugam(size(Heq{t},1)+1:end,:);
                for i = 1:N
                    dU{t,i} = u(start_m{i}+1:start_m{i}+m{i});
                    dL{t,i} = LAMLAM{t,i}*[1;dX{Lid{t}}];
                end
                dX{t+1} = F{t}*[1;dX{t};u];
            end
            mugam =  Ll{T+1}*[1;dX{Lid{T+1}}];    
            dM{T+1} = mugam(1:size(Heq{T+1},1),:);
            dG{T+1} = mugam(size(Heq{T+1},1)+1:end,:);
            dP = 0;
            flag = 1;
        end
    end
end