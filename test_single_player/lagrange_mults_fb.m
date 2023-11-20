n = 4;
m = 2;

Qa = rand(2*n,2*n);
Qa = Qa'*Qa + eye(2*n);
Qat = Qa(1:n,:);
Qad = Qa(n+1:end,:);
qa = rand(2*n,1);
qat = qa(1:n);
qad = qa(n+1:end);

Qb = rand(2*n,2*n);
Qb = Qb'*Qb + eye(2*n);
Qbt = Qb(1:n,:);
Qbd = Qb(n+1:end,:);
qb = rand(2*n,1);
qbt = qb(1:n);
qbd = qb(n+1:end);

Aa = rand(m,n);
Ab = rand(m,n);
ba = rand(m,1);
bb = rand(m,1);

Zmn = zeros(m,n);
Znm = zeros(n,m);
Zm = zeros(m,m); 




Mog = [Qat, Aa', Znm;
       Qbd, Znm, Ab';
       Aa,Zmn, Zm, Zm;
       Zmn, Ab,Zm, Zm];
mog = [qat; qbd; ba; bb];
soln = Mog\mog;

ta = soln(1:n);
tb = soln(n+1:2*n);
lama = soln(2*n+1:2*n+m);
lamb = soln(2*n+m+1:end);

cost_a1 = 0.5* [ta;tb]'*Qa*[ta;tb] + [ta;tb]'*qa;


Mi = inv(Mog)

dtadba = Mi(1:n,2*n+1:2*n+m);
dtadbb = Mi(1:n,2*n+m+1:end);

dtbdba = Mi(n+1:2*n,2*n+1:2*n+m);
dtbdbb = Mi(n+1:2*n,2*n+m+1:end);

dcadta = Qat * [ta;tb] - qat;
dcadtb = Qad * [ta;tb] - qad;

lamaa = -dtadba'*dcadta;
lamab = -dtadba'*dcadtb;

dcosta_dba = lamaa+lamab

ba2 = ba;
ba2(1) = ba2(1) + 0.01;
mog2 = [qat; qbd; ba2; bb];
soln2 = Mog\mog2;
ta2 = soln2(1:n);
tb2 = soln2(n+1:2*n);
cost_a2 = 0.5* [ta2;tb2]'*Qa*[ta2;tb2] + [ta2;tb2]'*qa;

dcost = (cost_a2-cost_a1) / 0.01

