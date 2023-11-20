%% LcpQR solver

%% Problem data
Q = eye(2);
q = zeros(2,1);

R = eye(1);
r = zeros(1,1);

S = zeros(1,2);

A = [1 0.1; 0 1];
B = [0.01/2; 0.1];
C = [0.01/2; 0.1];
d = zeros(2,1);

H = [1 0];
h = [0];

T = 10;

%% Solve

P = Q;
p = q;

G{T} = [];
g{T} = 

for t = T-1:-1:1
    
end

