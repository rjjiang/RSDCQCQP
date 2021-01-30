function result = eigqcqp(D1,D2,P2,b1,b2,ylb,yub,zlb,zub,L,d,m,n)%
model.obj = [b1;zeros(n,1)];
model.Q = sparse(blkdiag(D1,zeros(n)));
model.modelsense = 'min';

% nonconvex qc: x'*Q2*x +  c2'*x + d  <= 0
model.quadcon(1).Qc =  sparse(blkdiag(zeros(n),D2));
model.quadcon(1).q = [b2;zeros(n,1)];
model.quadcon(1).rhs = d;%result.x(n+1:2*n)'*D2*result.x(n+1:2*n)+b2'*result.x(1:n)+d
model.quadcon(1).sense = '<';
model.lb = [ylb;zlb];
model.ub = [yub;zub];

model.A = sparse([ [L zeros(m,n)]; [-eye(n) P2]]);
model.rhs = [ones(m,1); zeros(n,1)];
model.sense =  [repmat('<',m,1);repmat('=',n,1)];
% The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;
% if n==10
%    params.TimeLimit = '60';
% else
     params.TimeLimit = '300';    
% end
model.vtype = repmat('C',2*n,1);
result = gurobi(model, params);