function [D1, D2,newV, rcQ, flag]= RSDC(Q1,Q2,n)
flag = 0; %flag=1 if Q1 and Q2 are already SDC
Q12   = Q1\Q2;
[V, D] = eig(Q12);
d = diag(D);
complexQ = find(abs(imag(d))>1e-12);
realQ = setdiff((1:n)',complexQ);
if isempty(complexQ) %Q1 Q2 are already SD in this case
    VQV1 = V'*Q1*V; %return unnormalized
    VQV2 = V'*Q2*V;
    D1 = VQV1; %return unnormalized
    D2 = VQV2;
    rcQ = n;
    P1 = V;
    P2 = eye(n);
    flag = 1;
    fprintf("already SD\n");
    return;
elseif isempty(realQ) %Q1 Q2 are not SD in this case
    d = complexQ;
else
    d = d([realQ; complexQ]);
    V = V(:,[realQ; complexQ]);
end

D = diag(d);
%[Vnew,Dnew] = cdf2rdf(V,D);
%aug_Q1 =
%aug_Q2 =
rcQ = length(realQ);
% for i=1:rcQ
%     normalization = sqrt(abs(V(:, i)'*Q1*V(:, i)));
%     V(:, i) = V(:, i) /normalization;
% end
i = rcQ+1;
while (i<n+1)
    
    V(:,i:i+1) = [(V(:,i) + V(:,i+1)) / 2, 1i * (V(:,i) - V(:,i+1)) / 2];
    beta = V(:,i)'*Q1*V(:,i);
    alpha = V(:,i)'*Q1*V(:,i+1);
    if beta ~= 0
        y = beta / sqrt(2 * (alpha ^ 2 + beta ^ 2) * (alpha + sqrt(alpha^2 + beta^2)));
        x = (y / beta) * (-alpha - sqrt(alpha ^ 2 + beta ^ 2));
    elseif (beta == 0) && (alpha < 0)
        y = 0;
        x = sqrt(alpha);
    else
        x = 0;
        y = sqrt(-alpha);
    end
    V(:,i:i+1) = [x * V(:, i) + y * V(:, i + 1), -y * V(:, i) + x * V(:, i + 1)];
    i = i+2;
    
end
VQV1 = V'*Q1*V; %return unnormalized
VQV2 = V'*Q2*V;
if isempty(complexQ)~=1
    lcQ = length(complexQ);
    
    fghMatrix = (zeros(lcQ+1));
    rhs = zeros(lcQ+1,1);
    
    
    b = zeros(lcQ/2,1); %real part of eig
    a = zeros(lcQ/2,1); %image part of eig
    for i=1:(lcQ/2)
        a(i) = real(VQV2(n-lcQ+2*i-1,n-lcQ+2*i-1));
        b(i) = real(VQV2(n-lcQ+2*i-1,n-lcQ+2*i));
    end
    
    %dettt = 0;
    %while( abs(dettt) < 1e-12 | abs(dettt)>1e12)
    
    sigmas = ( linspace(-1,1,lcQ+1)');
    %sigmas = [a+b; a-b; 0];
    % sigmas = rand(lcQ+1,1)+0.5;
    i = 1;
    for j=1:lcQ+1
        fghMatrix(j,lcQ+1) =  (1) ; %%h
    end
    rhs = sigmas .*fghMatrix(:,lcQ+1);
    
    
    while (i <= lcQ/2)
        for j=1:lcQ+1
            fghMatrix(j,2*i-1) = ( 1 / ((b(i) - sigmas(j))^2 + a(i)^2 )); %f
            fghMatrix(j,2*i) = (sigmas(j) * fghMatrix(j,2*i-1));%g
        end
        i = i + 1;
    end
%%
    sol = fghMatrix\rhs;
    sol = real(sol);

    t = sol(1:2:lcQ-1);
    s = sol(2:2:lcQ);%2xy
    x = zeros(lcQ/2,1);
    y = zeros(lcQ/2,1);
    for i = 1:(lcQ/2)
        if t(i) == 0 %seems impossible (prob 0) since FGH matrix is invertible and rhs is nonzero
            if a(i) == 0
                fprintf('error: a(i)=0\n') ; %seems impossible since FGH matrix is invertible and rhs is nonzero
            elseif s(i)/a(i)>=0
                y(i) = sqrt(s(i)/a(i));
                x(i) = 0;
            else
                x(i) = sqrt(-s(i)/a(i));
                y(i) = 0;
            end
        else
            temp = (t(i) +  s(i) * b(i)) / a(i);
            y(i) = sqrt((temp + sqrt(temp^2 + s(i)^2 )) / 2);
            x(i) = s(i) / (2 * y(i));
        end
    end
    
    tcol = zeros(n,1);
    tcol(n-lcQ+1:2:n-1) = x;
    tcol(n-lcQ+2:2:n) = y;
    trow = tcol';
    newQ1 = [inv(V')*VQV1(1:n,1:n)*inv(V) zeros(n,1) ;zeros(n,1)'  1];
    newQ2 = [inv(V')*VQV2(1:n,1:n)*inv(V) inv(V')*tcol; trow*inv(V) sol(lcQ+1)];
    %resQ12 = inv(resQ1)*resQ2;
    newQ12 = newQ1\newQ2;
    [newV newD] = eig(newQ12);
    D1 = newV'*newQ1*newV;
    D2 = newV'*newQ2*newV;
%     for i=1:n+1
%         normalization = sqrt(abs(newV(:, i)'*resQ1*newV(:, i)));
%         newV(:, i) = newV(:, i) /normalization;
%     end
%    newVQV1 = newV'*newQ1*newV;
%    newVQV2 = newV'*newQ2*newV;
    fprintf('error=%e\n',norm(D1-diag(diag(D1)))+norm(D2-diag(diag(D2))));
end







