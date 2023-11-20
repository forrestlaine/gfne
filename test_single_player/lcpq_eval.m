%lcpq eval
clear; clc; close all;

dim_x = 2;
dim_u = 1;
dim_l = 1;

% Params for weird behavior
Q = diag([1,2]);
q = [0;10];

Qf = diag([1,2]);
qf = [0;-10];

H = [0 1];
h = 0;

% Q = diag([2,0]);
% q = [4;0];
% 
% Qf = diag([2,0]);
% qf = [-4;0];
% 
% H = [1 0];
% h = 2;


R = .1*eye(1);
r = [0];
RQ = [0. 0.];

eig([Q RQ'; RQ R])

T = 0.1;

A = [1 T; 0 1];
B = [T*T/2; T];
C = [T*T/2; T];
d = [0;0];


M = [0 (q+A'*(Qf*d+qf))' (r+B'*(Qf*d+qf))' (C'*(Qf*d+qf))';
    (q+A'*(Qf*d+qf)) Q+A'*Qf*A RQ'+A'*Qf*B A'*Qf*C;
    (r+B'*(Qf*d+qf)) RQ+B'*Qf*A R+B'*Qf*B B'*Qf*C;
    (C'*(Qf*d+qf)) C'*Qf*A C'*Qf*B C'*Qf*C];

Mshad = zeros(size(M));
Mshad(end,:) = M(end,:);
Mshad(:,end) = M(:,end);

Gx = H*A;
Gu = H*B;
g = H*d + h;

Lx = -(H*C)\(H*A);
Lu = -(H*C)\(H*B);
l = -(H*C)\(H*d+h);

p1 = [eye(1+dim_x+dim_u); 
      l Lx Lu];  
  
p2 = [eye(1+dim_x+dim_u);
      zeros(dim_l,1+dim_x+dim_u)];


M1 = p1'*M*p1;
M2 = p2'*M*p2;
Mshad1 = p1'*Mshad*p1;
Malt = p2'*Mshad*p2;


K1 = -M1(dim_x+2:end,dim_x+2:end)\M1(dim_x+2:end,1:dim_x+1);
K2 = -M2(dim_x+2:end,dim_x+2:end)\M2(dim_x+2:end,1:dim_x+1);

pk1 = [eye(1+dim_x);K1];
pk2 = [eye(1+dim_x);K2];
MM1 = pk1'*M1*pk1;
MM2 = pk2'*M2*pk2;

DM = MM1-MM2;

divider_a = [l Lx Lu]*pk1;
divider_b = [l Lx Lu]*pk2;

num_pts = 400;
X1 = linspace(-10,10,num_pts);
X2 = linspace(-20,20,num_pts);

check_out = false;
for i = 1:num_pts
    for j = 1:num_pts
        con = [1;X1(i);X2(j)];
        
        dividing_pt = -Gu\(Gx*[X1(i);X2(j)] + g);
        
        dividing_val = 0.5*[con;dividing_pt]'*M1*[con;dividing_pt];
        dividing_val2 = 0.5*[con;dividing_pt]'*M2*[con;dividing_pt];
        if abs(dividing_val-dividing_val2) > 1e-6
            disp('wtf')
        end
        if abs(Gx*[X1(i);X2(j)] + Gu*dividing_pt + g) > 1e-6
            disp('wtf2')
        end
        
        a1 = 0.5*M1(4,4);
        b1 = M1(4,1:3)*con;
        a2 = 0.5*M2(4,4);
        b2 = M2(4,1:3)*con;
        min1 = -b1/(2*a1);
        min2 = -b2/(2*a2);
        
        val_min_1 = 0.5*[con;min1]'*M1*[con;min1];
        val_min_2 = 0.5*[con;min2]'*M2*[con;min2];

        if min1 < dividing_pt
            left_min = min1;
            left_val = val_min_1;
        else
            left_min = dividing_pt;
            left_val = dividing_val;
        end
        if min2 > dividing_pt
            right_min = min2;
            right_val = val_min_2;
        else
            right_min = dividing_pt;
            right_val = dividing_val;
        end
        if left_val < right_val
            Z(i,j) = left_val;
            u_opt(i,j) = left_min;
        else
            Z(i,j) = right_val;
            u_opt(i,j) = right_min;
        end
        Z1(i,j) = val_min_1;
        Z2(i,j) = val_min_2;
        
        if Gx*[X1(i);X2(j)]+Gu*u_opt(i,j)+g < 0
%             if abs(u_opt(i,j)) > 0.1
%                 check_out = true;
%                 disp('why control here?')
%             end
            lam_opt(i,j) = Lx*[X1(i);X2(j)] + l + Lu*u_opt(i,j);
        else
            lam_opt(i,j) = 0;
        end
        
        if abs(X1(i)) < 0.25 && abs(X2(j)+0.5) < 0.25
            check_out = false;
        end
%         if abs(u_opt(i,j)) > 0.1
%             if abs(Lx*[X1(i);X2(j)] + l) > 0.1
%                 check_out = true;
%                 disp('what the actual fuck');
%             end
%         end
        x_next(i,j,1:2) = A*[X1(i);X2(j)] + B*u_opt(i,j) + C*lam_opt(i,j) + d;
        x_n = [x_next(i,j,1); x_next(i,j,2)];
        constraint_resid(i,j) = H*x_n + h;
        
        total_cost(i,j) = 0.5*(x_n'*Qf*x_n + u_opt(i,j)'*R*u_opt(i,j) + [X1(i);X2(j)]'*Q*[X1(i);X2(j)]);
        diff_cost(i,j) = 0.5*con'*(MM1-MM2)*con;
        if divider_a*con > 0 && divider_b*con > 0
            region(i,j) = 1;
        end
        if divider_a*con <= 0 && divider_b*con <= 0
            region(i,j) = 2;
        end
        if divider_a*con <= 0 && divider_b*con > 0
            region(i,j) = 3;
        end
        if divider_a*con > 0 && divider_b*con <= 0
            region(i,j) = 4;
%             if 0.5*(divider_a+divider_b)*con > 0
%                 u_opt(i,j) = K2*con;
%                 Z(i,j) = 0.5*con'*MM2*con;
%             else
%                 u_opt(i,j) = K2*con;
%                 Z(i,j) = 0.5*con'*MM2*con;
%             end
            
        end
        if (check_out)
            span = 2*abs(u_opt(i,j)-dividing_pt);
            if span < 10
                span = 10;
            end
            U1 = linspace(dividing_pt-span,dividing_pt,100);
            U2 = linspace(dividing_pt,dividing_pt+span,100);

            for k = 1:100
                y1(k) = 0.5*[con;U1(k)]'*M1*[con;U1(k)];
                y2(k) = 0.5*[con;U2(k)]'*M2*[con;U2(k)];   
                y3(k) = 0.5*[con;U1(k)]'*M2*[con;U1(k)];
                y4(k) = 0.5*[con;U2(k)]'*M1*[con;U2(k)];   
            end
            figure; hold on;
            plot(U1,y1);
            plot(U2,y2);
%             plot(U1,y3);
%             plot(U2,y4);
            plot(u_opt(i,j),Z(i,j),'or', 'MarkerSize', 4, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
            plot(min1,val_min_1,'or','MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');
            plot(min2,val_min_2,'or','MarkerSize', 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'b');   
            X1(i)
            X2(j)
            close all;
        end
    end
end
%%
figure;
contour(X1,X2,Z',100);
hold on;
% figure;
% contour(X1,X2,u_opt',100);
% figure;
% hold on;
% figure;
fimplicit(@(x1,x2) [1; x1;x2]'*(MM1-MM2)*[1;x1;x2],'-r','linewidth',4)


[V,D] = eig(DM);
% Dalt = D;
% Dalt(2,2) = 1;
% Dalt(1,1) = -1;
% DMalt = V*Dalt*V';
% fimplicit(@(x1,x2) [1; x1;x2]'*DMalt*[1;x1;x2],'-g','linewidth',2)
figure; 
contour(X1,X2,diff_cost',100);
hold on;
fimplicit(@(x1,x2) [1; x1;x2]'*(MM1-MM2)*[1;x1;x2],'-r','linewidth',2)
% lam = diag(D);
% utotheta = [1 0; 0 1; sqrt(-D(1,1)/D(3,3)),0];
% coefs = V(1,:)*utotheta;
% utotheta2 = [1 0; 0 1; -sqrt(-D(1,1)/D(3,3)),0];
% coefs2 = V(1,:)*utotheta2;
% slope2 = -coefs2(1)/coefs2(2);
% bias2 = 1/coefs2(2);
% 
% U1 = linspace(-3,3,100);
% U2 = slope*U1+bias;
% U22 = slope2*U1+bias2;
% lineX = V*utotheta*[U1;U2];
% lineX2 = V*utotheta2*[U1;U22];
% plot(lineX(2,:),lineX(3,:),'-g','linewidth',2);
% plot(lineX2(2,:),lineX2(3,:),'-g','linewidth',2);
% figure;
% fill(X1,X2,region1');
% figure;
% contour(X1,X2,region',100);
% figure;
% contour(X1,X2,u_opt',100);
% fimplicit(@(x1,x2) (MM1(2,:)-MM2(2,:))*[1;x1;x2],'-k','linewidth',4)
% fimplicit(@(x1,x2) (MM1(3,:)-MM2(3,:))*[1;x1;x2],'-k','linewidth',4)
% fimplicit(@(x1,x2) divider_a*[1; x1;x2],'-g','linewidth',4)
% fimplicit(@(x1,x2) divider_b*[1; x1;x2],'-b','linewidth',4)
% fimplicit(@(x1,x2) 0.5*(divider_a+divider_b)*[1; x1;x2],'-k','linewidth',4)
% figure;
% hold on;
% contour(X1,X2,Z1',200);
% figure;
% hold on;
% contour(X1,X2,Z2',200);


% figure;
% contour(X1,X2,u_opt',400);
% figure; 
% contour(X1,X2,lam_opt',400);
% figure;
% contour(X1,X2,x_next(:,:,1)',400);

% figure;
% contour(X1,X2,constraint_resid',400);
 
% Q1 = [3 1; 1 2];
% q1 = [2;2];
% 
% Q2 = [2 1; 1 3];
% q2 = [-2;2];
% 
% X = linspace(-30,30,200);
% Y = linspace(-30,30,200);
% 
% for i = 1:200
%     for j = 1:200
%         pt = [X(i);Y(j)];
%         Z1(i,j) = 0.5*pt'*Q1*pt + pt'*q1;
%         Z2(i,j) = 0.5*pt'*Q2*pt + pt'*q2;
%     end
% end
% Z3 = min(Z1,Z2);
% contour(X,Y,Z3');
% hold on;
% fimplicit(@(x1,x2) 0.5*[x1;x2]'*(Q1-Q2)*[x1;x2] + [x1;x2]'*(q1-q2),'-b')
% % hold on;
% % contour(X,Y,Z2)