function result = oriqcqp(Q1,Q2,c1,c2,lb,ub,d,A,b,n)
% Objective function min x'*Q1*x +  c1'*x 
model.obj = c1;
model.Q = sparse(Q1);
model.modelsense = 'min';

% nonconvex qc: x'*Q2*x +  c2'*x + d  <= 0
model.quadcon(1).Qc =  sparse(Q2);
model.quadcon(1).q = c2;
model.quadcon(1).rhs = d;
model.quadcon(1).sense = '<';
model.lb = lb;
model.ub = ub;



%linear constraints: Ax<=b
model.A = sparse(A);
model.rhs = b;


% The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;

%if n==10
 %   params.TimeLimit = '60';
%else
   params.TimeLimit = '300';    
%end
model.vtype = repmat('C',n,1);
result = gurobi(model, params);