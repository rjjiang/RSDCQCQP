%main
%solve qcqp that the quadratic forms are SDC
% solve    min  x'*Q1*x+c1'*x
%          s.t. x'*Q2*x+c2'*x + d<=0
%               x'*x<=1
%               A*x<=b
clear;clc

for n = [10 15 20]
    fid = fopen(['sdqcqp_',num2str(n),'.txt'],'wt');
    fprintf(fid,'numerical comparison\n')
    fprintf(fid,'------------------------------------------------------------------\n')
    i=1;
    while (i<=5)
 %       try
            m=100;
            A = randn( n);  V = orth( A.' ).';
            d1 = diag(sign(randn(n,1)));
            d2 = diag(randn(n,1));
            Q1 = V'*d1*V;
            Q1 = (Q1+Q1')/2;
            Q2 = V'*d2*V;
            Q2 = (Q2+Q2')/2;
            b1 = 2* randn(n,1);
            b2 = 2* randn(n,1);
            x0 = zeros(n,1);
            d = 1;
            L = randn(m,n);
            b = [L*x0 + 1];
            
            Q12 =  Q1\Q2;
            [P1,D] = eig(Q12);
            D1 = P1'*Q1*P1;
            D2 = P1'*Q2*P1;
            %% oriQCQP
            [xlb,xub] = bound(L,b,m,n);
            resultori = oriqcqp(Q1,Q2,b1,b2,xlb,xub,d,L,b,n);
            
            
            %% 1RSDC
            if norm(abs(D1-diag(diag(D1))))+norm(abs(D2-diag(diag(D2)))) <1e-6
                fprintf('SDC is done\n');
            else
                fprintf('not SDC\n');continue
            end
            hc1 = P1' * b1;
            hc2 = P1' * b2;
            LP1 = L*P1;
            sD1 = sparse(diag(diag(D1)));
            sD2 = sparse(diag(diag(D2)));
            [wlb,wub] = bound(LP1,b,m,n);
            resultsdc =  oriqcqp(sD1,sD2,hc1,hc2,wlb,wub,d,LP1,b,n);% DCqpGu(sD1,sD2,hc1,hc2,wlb,wub,d,tildeAP2,b,P2np1,m,n+1);% with one linear equation
            
            %% eigQCQP
            [V1, eD1] = eig(Q1);%   V1'*Q1*V1 = eD1  V1y=x
            [V2, eD2] = eig(V1'*Q2*V1);
            [ylb, yub] = boundyz(V1,L,b,m,n);
            [zlb, zub] = boundyz(V1*V2,L,b,m,n);
            LV1 = L*V1;
            bv1 = V1'*b1;
            bv2 = V1'*b2;
            resulteig = eigqcqp(eD1,eD2,V2,bv1,bv2,ylb,yub,zlb,zub,LV1,d,m,n);% with one linear equation   x=V1*resulteig.x(1:n)
            %%
            fprintf(fid,'i=%d\t\n',i);
            fprintf(fid,'\n oriQCQP:\t &%1.4e\t &%5.2f\n',resultori.objval,resultori.runtime);
            fprintf(fid,'sdcQCQP:\t &%1.4e\t &%5.2f\t &%1.4e\n',resultsdc.objval,resultsdc.runtime,norm(P1)*norm(P1));
            fprintf(fid,'eigQCQP:\t &%1.4e\t &%5.2f\n',resulteig.objval,resulteig.runtime);
            fprintf(fid,'------------------------------------------------------------------\n');
            i=i+1;
%        catch
%            i=i+1;
%        end
    end
    fclose(fid);
end