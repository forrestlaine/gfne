function [z, retcode] = LMCP(M,q,l,u,z0)
% LMCP Solves linear mixed complementarity problems:
%   Find z such that
%   l_i < z_i < u_i --> (Mz+q)_i  = 0
%   l_i = z_i       --> (Mz+q)_i >= 0
%         z_i = u_i --> (Mz+q)_i <= 0
% The algorithm is based off of what is described in 
% Path: A nonmonotone ...
% Optional: z0, the initial guess of solution. 
% If z0 not provided, it is initialized as follows:
%   -if l_i is finite, z0_i = l_i
%   -else if u_i is finite, z0_i = u_i
%   -else, z0_i = 0
n = numel(q);
if (n~=numel(l) ||  n~=numel(u) || n~=size(M,1) || n~=size(M,2))
    retcode = -1; % data dimension mismatch
    z = NaN*ones(n,1);
elseif any(u < l)
    retcode = -2; % infeasible problem
    z= NaN*ones(n,1);
else
    if nargin < 5
        unbounded_below = isinf(l);
        unbounded_above = isinf(u);
        z0 = zeros(n,1);
        for i = 1:n
            if unbounded_below(i)
                if unbounded_above(i)
                    z0(i) = 0;
                else
                    z0(i) = u(i);
                end
            else
                z0(i) = l(i);
            end
        end   
        x0 = z0;
    else
        x0 = z0;
        z0 = max(min(x0,u),l);
    end
            
    r = M*z0+q + x0-z0;
    
    A = [M -eye(n) eye(n)  r];
    a = [-q + r];
    m = size(A,2);
    
    aa = a;
    basis_set = [];
    for i = 1:n
        if (z0(i)-l(i) <= 1e-6)
            basis_set = [basis_set, n+i];
            aa = aa - l(i)*M(:,i);
        elseif (u(i)-z0(i) <= 1e-6)
            basis_set = [basis_set, n+n+i];
            aa = aa - u(i)*M(:,i);
        else
            basis_set = [basis_set, i];
        end
    end
    basis_set = sort(basis_set);
    
    B = A(:,basis_set);
    b = B\aa;

    entering_index = m;
    entering_value = 0;
    opposing_value = 1;
    dir = sign(opposing_value-entering_value);
    
    for iters = 1:100
        ev = B\A(:,entering_index);
        ratios = zeros(n,1);
        ut = zeros(n,1);
        lt = zeros(n,1);
        for i = 1:n
            if ev(i) == 0
                ratios(i) = Inf;
            elseif basis_set(i) <= n
                ut(i) = (b(i)-u(i)-ev(i)*entering_value)/(ev(i)*dir);
                lt(i) = (b(i)-l(i)-ev(i)*entering_value)/(ev(i)*dir);
                ratios(i) = max(ut(i),lt(i));
            elseif basis_set(i) <= 3*n
                ut(i) = (b(i)-inf-ev(i)*entering_value)/(ev(i)*dir);
                lt(i) = (b(i)-0-ev(i)*entering_value)/(ev(i)*dir);
                ratios(i) = max(ut(i),lt(i));
            else
                ut(i) = (b(i)-1-ev(i)*entering_value)/(ev(i)*dir);
                lt(i) = (b(i)-(-inf)-ev(i)*entering_value)/(ev(i)*dir);
                if ut(i) > 0
                    ratios(i) = ut(i);
                else
                    ratios(i) = inf;
                end
%                 ratios(i) = max(ut(i),lt(i));
            end
        end
        [min_t, min_index] = min(ratios, [], 1);
        if (isinf(min_t) && iters > 1)
            retcode = -3;
            z = NaN*ones(n,1);
            break;
        end
        if (min_t >= (opposing_value-entering_value)/dir)
            exiting_index = entering_index;
            exiting_down = -dir;
        else
            
            exiting_index = basis_set(min_index);
            exiting_down = lt(min_index) > ut(min_index);
      
        end
        basis_set = [basis_set, entering_index];
        basis_set = sort(basis_set(basis_set~=exiting_index));
        B = A(:,basis_set);
        aa = a;
        for j = basis_set
            if j > n && j <= 2*n
                aa = aa - l(j-n)*M(:,j-n);
            elseif j > 2*n && j <= 3*n
                aa = aa - u(j-2*n)*M(:,j-2*n);
            end
        end
        if ~ismember(basis_set,m)
            aa = aa - r;
        end
        
        b = B\aa;

        
        if exiting_index == m
            zwv = zeros(3*n,1);
            zwv(basis_set) = b;
            for i = 1:n
                if ismember(i+n, basis_set)
                    zwv(i) = l(i);
                elseif ismember(i+2*n, basis_set)
                    zwv(i) = u(i);
                end
            end            
            z = zwv(1:n);
            retcode = 1;
            break;
        else
            if exiting_index <= n
                entering_value = 0;
                opposing_value = inf;
                if exiting_down
                    entering_index = exiting_index+n;
                else
                    entering_index = exiting_index+2*n;
                end
            elseif exiting_index <= 2*n
                entering_index = exiting_index - n;
                entering_value = l(entering_index);
                opposing_value = u(entering_index);
            else 
                entering_index = exiting_index - 2*n;
                entering_value = l(entering_index);
                opposing_value = u(entering_index);
            end
        end

    end
end

