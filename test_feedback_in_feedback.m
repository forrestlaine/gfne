% Test feedback under feedback

Qa = randn(3,3);
Ra = randn(2,2);
Ra = Ra'*Ra+eye(2);
Qa = Qa'*Qa;
Sa = 0.1*randn(3,2);

Qb = randn(3,3);
Rb = randn(2,2);
Rb = Rb'*Rb+eye(2);
Qb = Qb'*Qb;
Sb = 0.1*randn(3,2);

Pa = Qa;
Pb = Qb;

A = randn(3,3) + eye(3);
Ba = randn(3,1);
Bb = randn(3,1);

Ha = [Qa Sa;Sa' Ra] + [A Ba Bb]'*Pa*[A Ba Bb];
Hb = [Qb Sb;Sb' Rb] + [A Ba Bb]'*Pb*[A Ba Bb];

Mopt = [Ha(4,:); Hb(5,:)];
Kopt = -Mopt(:,4:5)\Mopt(:,1:3);

Krand = randn(2,3);

pol = [eye(3) zeros(3,2); Krand eye(2)];
Har = pol'*Ha*pol;
Hbr = pol'*Hb*pol;

Mrand = [Har(4,:); Hbr(5,:)];
Krandopt = -Mrand(:,4:5)\Mrand(:,1:3);
