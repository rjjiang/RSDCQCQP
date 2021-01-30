function [ylb,yub] = boundw(A,b,a,m,n)
[x,y]=size(a);%
ylb = -inf*ones(n,1);
yub = inf*ones(n,1);
for i = 1:n
    c = zeros(n,1);
    c(i) = 1;
    model.obj = c;
    model.modelsense = 'Min';
    
    %linear constraints: Ax<=b
    % model.A = sparse([A; ae']);
    % model.rhs = [b; 0];
    % model.sense =  [repmat('<',m,1); '='];
    model.A = sparse([A;a']);
    model.rhs = [b; zeros(y,1)];%improve stability
    
    model.sense =  [repmat('<',m,1);repmat('=',y,1)];
    % model.sense =  [repmat('<',m+2,1)];
    model.lb = -inf*ones(n,1);
    model.ub = inf*ones(n,1);
    % The problem is non-convex,
    % we need to set the parameter 'NonConvex' in order to solve it.
    %params.NonConvex = 2;
    params.TimeLimit = '60';params.OutputFlag = 0;
    model.vtype = repmat('C',n,1);
    result = gurobi(model, params);
    if strcmp(result.status,'OPTIMAL') ==1
        ylb(i) = result.objval;
    end
    
    model.modelsense = 'max';
    result = gurobi(model, params);
    if strcmp(result.status,'OPTIMAL') ==1
        yub(i) = result.objval;
    end
end

end