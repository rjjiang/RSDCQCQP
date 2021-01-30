function [D1, D2, P, lrQ, flag]= RSDC2(Q1,Q2,n)
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
    lrQ = n;
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
lrQ = length(realQ);
% for i=1:lrQ
%     normalization = sqrt(abs(V(:, i)'*Q1*V(:, i)));
%     V(:, i) = V(:, i) /normalization;
% end
i = lrQ+1;
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
    
    fghMatrix = (zeros(lcQ/2+1));
    rhs = zeros(lcQ/2+1,1);
    
    
    reQ2 = zeros(lcQ/2,1); %real part of eig
    imQ2 = zeros(lcQ/2,1); %image part of eig
    lambda = zeros(lcQ/2,1);
    for i=1:(lcQ/2)
        imQ2(i) = -real(VQV2(n-lcQ+2*i-1,n-lcQ+2*i-1));
        reQ2(i) = real(VQV2(n-lcQ+2*i-1,n-lcQ+2*i));
        lambda(i) = reQ2(i)+1i*imQ2(i);
    end
    
    sigmas = ( linspace(-1,1,lcQ/2+1)');
    
    i = 1;
    for j=1:lcQ/2+1
        fghMatrix(j,lcQ/2+1) =  (1) ; %%h
    end
    rhs = sigmas .*fghMatrix(:,lcQ/2+1);
    while (i <= lcQ/2)
        for j=1:lcQ/2+1
            fghMatrix(j,i) = -( 1 / (lambda(i)-sigmas(j)));
        end
        i = i + 1;
    end
    sol = fghMatrix\rhs;
    for i=1:lcQ/2
        sol(i) = sqrt(sol(i));
    end
%% we shold have two same eigenvalues...
    
    
    alam = real(sol);
    blam = imag(sol);
    
    tcol = zeros(lrQ,2);
    for i=1:lcQ/2
        tcol = [tcol; [-blam(i) alam(i);alam(i) blam(i)]];
    end
    trow = tcol';  
    resQ1 = [inv(V')*VQV1*inv(V) zeros(n,2) ;zeros(n,2)'  [0 1;1 0]];
    resQ2 = [inv(V')*VQV2*inv(V) inv(V')*tcol; trow*inv(V) [-imag(sol(lcQ/2+1)) real(sol(lcQ/2+1)) ;real(sol(lcQ/2+1)) imag(sol(lcQ/2+1))]];
    %resQ12 = inv(resQ1)*resQ2;
    resQ12 = resQ1\resQ2;
    [resV resD] = eig(resQ12);
    RresV = resV;
    vnorm = zeros(n+2,1);
    for i=1:n+2
        vnorm(i) = norm(imag(resV(:,i)));
    end
    complexv = find(abs(vnorm)>1e-8);
    realv = setdiff((1:n+2)',complexv);
    RresV = real(resV(:,realv));
    temp = resV(:,complexv);
    for i = 1:2:length(complexv)-1
            RresV = [RresV, real(temp(:,i)+temp(:,i+1))];
            RresV = [RresV,imag(temp(:,i)-temp(:,i+1))];
    end
    
   
    resVQV1 = RresV'*resQ1*RresV;
    resVQV2 = RresV'*resQ2*RresV;
    [eV1 eD]=eig( (resVQV1'+resVQV1)/2);%in case that multiplicity of real eigenvalues is >=2
    RresV= RresV*eV1; 
    
%     for i=1:lcQ+2
%         normalization = sqrt(abs(ResV(:, i)'*resQ1*ResV(:, i)));
%         ResV(:, i) = ResV(:, i) /normalization;
%     end
    resVQV1 = RresV'*resQ1*RresV;
    resVQV2 = RresV'*resQ2*RresV;
    [eV2 eD]=eig( (resVQV2'+resVQV2)/2);%in case that multiplicity of real eigenvalues is >=2
    RresV= RresV*eV2; 
    D1 = RresV'*resQ1*RresV;
    D2 = RresV'*resQ2*RresV;
    P = RresV;
%     D1 = blkdiag(VQV1(1:n-lcQ,1:n-lcQ), resVQV1 );
%     D2 = blkdiag(VQV2(1:n-lcQ,1:n-lcQ), resVQV2 );
%     P1 = V;
%     P2 = blkdiag(eye(n-lcQ),RresV);
%     p = [zeros(n-lcQ,2);tcol];
  fprintf('error=%e\n',norm(D1-diag(diag(D1)))+norm(D2-diag(diag(D2))));
end
