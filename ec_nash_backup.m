%% ec_nash_backup

function [K, P, G] = ec_nash_backup(F,G,P,H,Q,N)
    n = size(F,1); % shared state dim
    
    for i = 1:N
        m{i} = size(H,2) - n - 1; % private control dim
        HH = [H; 
        
        % G{i}*[1; x; ui] = 0;
        Gcx = G{i}(:,1:n+1);
        Gu = G{i}(:,2+n:end);
        [U,S,V] = svd(Gu);
        rk = rank(S,1e-6);
        U1 = U(:,1:rk);
        U2 = U(:,rk+1:end);
        S = S(1:rk,1:rk);
        V1 = V(:,1:rk);
        V2 = V(:,rk+1:end);
        
        
    end

end