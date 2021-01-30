function [ylb, yub] = boundyz(V,A,b,m,n)
ylb = zeros(n,1);
yub = zeros(n,1);
for i = 1:n
    c = zeros(2*n,1);
    c(n+i) = 1;
    model.obj = c;
    model.modelsense = 'Min';
    
    %linear constraints: Ax<=b
    % model.A = sparse([A; ae']);
    % model.rhs = [b; 0];
    % model.sense =  [repmat('<',m,1); '='];
     model.A = sparse([[A  zeros(m,n)];[-eye(n) V];[eye(n) -V]]);
     model.rhs = [b; 1e-8*ones(n,1);1e-8*ones(n,1)];%improve stability
     model.sense =  [repmat('<',m+2*n,1)];
%    model.A = sparse([[A  zeros(m,n)];[zeros(n) V]]);
%    model.rhs = [b; 0*ones(n,1)];%improve stability
%    model.sense =  [repmat('<',m,1);repmat('=',n,1)];
    
    model.lb = -inf*ones(2*n,1);
    model.ub = inf*ones(2*n,1);
    % The problem is non-convex,
    % we need to set the parameter 'NonConvex' in order to solve it.
    %params.NonConvex = 2;
    params.TimeLimit = '60';params.OutputFlag = 0;
    model.vtype = repmat('C',2*n,1);
    result = gurobi(model, params);
    ylb(i) = result.objval;
    
    model.modelsense = 'max';
    result = gurobi(model, params);
    yub(i) = result.objval;
end
    
end