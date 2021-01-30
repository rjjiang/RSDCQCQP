function result = twoqcqp(Q1,Q2,c1,c2,lb,ub,d,A,b,ae,m,n)
[~,y]=size(ae);%
model.obj = c1;
model.Q = sparse(Q1);
model.modelsense = 'min';

% nonconvex qc: x'*Q2*x +  c2'*x + d  <= 0
model.quadcon(1).Qc =  sparse(Q2);
model.quadcon(1).q = c2;
model.quadcon(1).rhs = d;%result.x'* Q2*result.x + c2'*result.x
model.quadcon(1).sense = '<';
model.lb = lb;
model.ub = ub;
% for i=1:m
%    b(i) = b(i) / norm(A(i,:));
%    A(i,:)  = A(i,:) / norm(A(i,:));
% end
%ae = ae / norm(ae);
model.A = sparse([A; ae']);
model.rhs = [b; zeros(y,1)];
model.sense =  [repmat('<',m,1); repmat('=',y,1)];
%model.A = sparse([A; ae';-ae']);
%model.rhs = [b; 1e-8;1e-8];
%model.sense =  [repmat('<',m+2,1)];
% The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;
%if n==11 || n==12
%   params.TimeLimit = '60';
%else
    params.TimeLimit = '300';    
%end
model.vtype = repmat('C',n,1);
result = gurobi(model, params);