function [ylb,yub] = bound(pA,b,m,n)
ylb = zeros(n,1);
yub = zeros(n,1);
for i = 1:n
    c = zeros(n,1);
    c(i) = 1;
    model.obj = c;
    model.modelsense = 'Min';
    
    %linear constraints: Ax<=b
    % model.A = sparse([A; ae']);
    % model.rhs = [b; 0];
    % model.sense =  [repmat('<',m,1); '='];
    model.A = sparse(pA);
    model.rhs = b;
    model.sense =  [repmat('<',m,1)];
    model.lb = -inf*ones(n,1);
    model.ub = inf*ones(n,1);
    % The problem is non-convex,
    % we need to set the parameter 'NonConvex' in order to solve it.
    %params.NonConvex = 2;
    params.TimeLimit = '60';params.OutputFlag = 0;
    model.vtype = repmat('C',n,1);
    result = gurobi(model, params);
    ylb(i) = result.objval;
    
    model.modelsense = 'max';
    result = gurobi(model, params);
    yub(i) = result.objval;
end
    
end