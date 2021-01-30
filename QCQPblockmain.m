%main
%solve qcqp that the quadratic forms are not SDC
% solve    min  x'*Q1*x+c1'*x
%          s.t. x'*Q2*x+c2'*x + d<=0
%               x'*x<=1
%               A*x<=b 
% Alex L. Wang and Rujun Jiang. New notions of simultaneous diagonalizability of quadratic forms with applications to QCQPs. arXiv preprint arXiv:2101.12141, 2021. code

clear;clc


for n = [10]
    for r =  [n-2]
        fid = fopen(['blockqcqp_n=',num2str(n),'_r=_',num2str(r),'.txt'],'wt');
        fprintf(fid,'numerical comparison\n');
        fprintf('------------------------------------------------------------------\n');
        fprintf(fid,'name\t objval\t runtime\t cond \t negative eigenvalues\n');
        i = 1;
        while (i <= 5)
     %       try
                m=100;
                A=randn( n);  V=orth( A.' ).';
                F = [0 1;1 0];
                d1 = [];
                d2 = [];
                if r>0
                    d1 = diag(sign(randn(r,1)));
                    d2 =  diag(randn(r,1));
                end
                for jj=1:(n-r)/2
                    d1 = blkdiag(d1,F);
                    a = randn(1);
                    b = randn(1);
                    FG = [b a;a -b];
                    d2 = blkdiag(d2,FG);
                end
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
                % x0 is a feasible solution
                start = tic;
                
 %% ori
                [xlb,xub] = bound(L,b,m,n);
                if isinf(sum(abs(xlb))+sum(abs(xub))) ==1
                    continue;
                end
                resultori = oriqcqp(Q1,Q2,b1,b2,xlb,xub,d,L,b,n);
                

%% 1RSDC 
                start = tic;
                [D1, D2,P,r, flag ]= RSDC(Q1,Q2,n);
                timesdc = toc(start);
                if flag == 0
                    hb1 = P' *  [ b1; 0];
                    hb2 = P' *  [b2; 0];
                    LP = [L,zeros(m,1)]*P;
                    Pnp1 = P(n+1,:)';%make row vector to column vector
                    sD1 = sparse(diag(diag(D1)));
                    sD2 = sparse(diag(diag(D2)));
                    [wlb,wub] = boundw(LP,b,Pnp1,m,n+1); % y = P2*resultsdc.x; x=P1*y(1:n)
                    result1rsdc = twoqcqp(sD1,sD2,hb1,hb2,wlb,wub,d,LP,b,Pnp1,m,n+1);% with one linear equation        
 %% 2RSDC
                start = tic;
                [D1, D2, Pr,r, flag ]= RSDC2(Q1,Q2,n);%temp=inv(Pr2')*D2*inv(Pr2);temp(1:n,1:n)-VQV2
                %VQV1 - Pr1'* Q1*Pr1
                timesdc = toc(start);
                neg2 = 0;
                for ii=1:n+2
                    if D1(ii,ii)<=0
                        neg2=neg2+1;
                    end
                    if D2(ii,ii)<=0
                        neg2=neg2+1;
                    end
                end
                if flag == 0
                    hb1 = Pr' *  [b1; 0;0];
                    hb2 = Pr' *  [b2; 0;0];
                    LP = [L,zeros(m,2)]*Pr;
                    Pnp1 = Pr(n+1:n+2,:)';%make row vector to column vector
                    sD1 = sparse(diag(diag(D1)));
                    sD2 = sparse(diag(diag(D2)));
                    [wlb,wub] = boundw(LP,b,Pnp1,m,n+2);
                    result2rsdc = twoqcqp(D1,D2,hb1,hb2,wlb,wub,d,LP,b,Pnp1,m,n+2);
                end
%% eigQCQP
                    [V1, eD1] = eig(Q1);%   V1'*Q1*V1 = eD1  V1y=x
                    [V2, eD2] = eig(V1'*Q2*V1);
                    [ylb, yub] = boundyz(V1,L,b,m,n);
                    [zlb, zub] = boundyz(V1*V2,L,b,m,n);
                    LV1 = L*V1;
                    bv1 = V1'*b1;
                    bv2 = V1'*b2;
                    resulteig = eigqcqp(eD1,eD2,V2,bv1,bv2,ylb,yub,zlb,zub,LV1,d,m,n);% with one linear equation   
%%
                    fprintf(fid,'i=%d\t real eigenvalue #=%d\n',i,r);
                    fprintf(fid,'(%d,%d) \t %1.4e\t &%5.2f\t& \n',n,n-r,resultori.objval,resultori.runtime);
                    fprintf(fid,' \t%1.3e\t &%5.2f\t& %1.2e&\n',result1rsdc.objval,result1rsdc.runtime, ...
                        norm(P)*norm(inv(P)));
                    fprintf(fid,' \t%1.3e\t &%5.2f\t& %1.2e&\n',result2rsdc.objval,result2rsdc.runtime, ...
                        norm(Pr)*norm(inv(Pr)));
                    fprintf(fid,'\t%1.3e\t &%5.2f\n',resulteig.objval,resulteig.runtime);
                    fprintf(fid,'------------------------------------------------------------------\n');
                    i=i+1;
                end
  %          catch
  %              i=i+1;
  %          end
        end
        fclose(fid);
    end
end