% test nested opt/ eq

n = 5;
m1 = 2;
m2 = 2;
nm = n+m1+m2;

F = randn(n,nm);

Q1 = randn(nm,nm);
Q1 = Q1'*Q1 + eye(nm);

Q2 = randn(nm,nm);
Q2 = Q2'*Q2 + eye(nm);


p = 3;
G = randn(p,n+m1+m2);


M1 = [Q1(1:n+m1,:), F(:,1:n+m1)', zeros(n+m1,n), G(:,1:n+m1)';
     Q2([1:n, n+m1+1:n+m1+m2],:), zeros(n+m2,n), F(:,[1:n, n+m1+1:n+m1+m2])', G(:,[1:n, n+m1+1:n+m1+m2])';
     F, zeros(n,n+n+p);
     G, zeros(p,n+n+p)];
 
M2 = [Q1(1:n+m1,:), F(:,1:n+m1)', zeros(n+m1,n), G(:,1:n+m1)' zeros(n+m1,p);
     Q2([1:n, n+m1+1:n+m1+m2],:), zeros(n+m2,n), F(:,[1:n, n+m1+1:n+m1+m2])', zeros(n+m2,p), G(:,[1:n, n+m1+1:n+m1+m2])';
     F, zeros(n,n+n+2*p);
     G, zeros(p,n+n+2*p);
     zeros(p,n+n+n+m1+m2), eye(p), -eye(p)];
 
b1 = randn(size(M1,1),1);
b2 = [b1; randn(p,1)];

d1 = M1\b1;
d2 = M2\b2;

