function [A, AL, S, Z, phi, alpha] = BPFA_SR(X,XL,K,iter, sampleA, A, phi, alpha)

[P, N] = size(X);
[PL, ~] = size(XL);
AL = [];
phiL = 1;

%hyper-parameters setting
%hyper-parameters for pi (usaully set b0 to a large number)
a0=1;
b0=N/8;
%hyper-parameters for phi
c0=1e-6;
d0=1e-6;
%hyper-parameters for alpha
e0=1e-6;
f0=1e-6;

%initialize the latent variables
if nargin<5, sampleA = true; end
if sampleA, A=randn(P,K); AL=randn(PL,K); end
if nargin<7, phi=1; end
if nargin<8, alpha=ones(1,K); end
S = randn(N,K);
Z = false(N,K);
Pi = ones(K,1);

Index = ones(size(X));
if sampleA, IndexL = ones(size(XL)); end

%Gibbs sampling part
for it = 1:iter
    X_k = Index.*X-Index.*(A*(Z.*S)');
    
    if sampleA
        X_kL = IndexL.*XL-IndexL.*(AL*(Z.*S)');
        
        %sample A
        for j=1:K
            X_k(:,Z(:,j)) = X_k(:,Z(:,j)) + sparse_mult(Index(:,Z(:,j)),A(:,j),S(Z(:,j),j));
            sig_A = 1./(P + phi*Index(:,Z(:,j))*(S(Z(:,j),j).^2));
            mu_A=phi*sig_A.*(X_k(:,Z(:,j))*S(Z(:,j),j));
            A(:,j)=sqrt(sig_A).*randn(P,1)+mu_A;
            X_k(:,Z(:,j)) = X_k(:,Z(:,j)) - sparse_mult(Index(:,Z(:,j)),A(:,j),S(Z(:,j),j));
        end
        info(it,X,A,S',Z',phi,alpha,'DH');

        for i=1:30
            for j=1:K
                X_kL(:,Z(:,j)) = X_kL(:,Z(:,j)) + sparse_mult(IndexL(:,Z(:,j)),AL(:,j),S(Z(:,j),j));
                sig_A = 1./(PL + phiL*IndexL(:,Z(:,j))*(S(Z(:,j),j).^2));
                mu_A=phiL*sig_A.*(X_kL(:,Z(:,j))*S(Z(:,j),j));
                AL(:,j)=sqrt(sig_A).*randn(PL,1)+mu_A;
                X_kL(:,Z(:,j)) = X_kL(:,Z(:,j)) - sparse_mult(IndexL(:,Z(:,j)),AL(:,j),S(Z(:,j),j));
            end
            ci=c0+1/2*sum(sum(Index));
            di=d0+1/2*sum(sum((IndexL.*XL-IndexL.*(AL*(Z.*S)')).^2));
            phiL=gamrnd(ci,1./di);
        end
        info(it,XL,AL,S',Z',phiL,alpha,'DL');
    end
    
    % sample Z
    for j=1:K
        X_k(:,Z(:,j)) = X_k(:,Z(:,j)) + sparse_mult(Index(:,Z(:,j)),A(:,j),S(Z(:,j),j));
        tempz1=((A(:,j).^2)'*Index)';
        tempz2=(A(:,j)'*X_k)';
        tempz=-(S(:,j).^2.*tempz1/2-S(:,j).*tempz2)*phi;
        tempz=exp(tempz)*Pi(j);
        Z(:,j)=rand(N,1)>((1-Pi(j))./(tempz+1-Pi(j)));
        X_k(:,Z(:,j)) = X_k(:,Z(:,j)) - sparse_mult(Index(:,Z(:,j)),A(:,j),S(Z(:,j),j));
    end
    info(it,X,A,S',Z',phi,alpha,'Z');
    
    % sample S
    for j=1:K
        X_k(:,Z(:,j)) = X_k(:,Z(:,j)) + sparse_mult(Index(:,Z(:,j)),A(:,j),S(Z(:,j),j));
        numz=nnz(Z(:,j));
        temps1=((A(:,j).^2)'*Index(:,Z(:,j)))';
        temps2=(A(:,j)'*X_k(:,Z(:,j)))';
        sig_s=1./(alpha(j)*ones(numz,1)+phi*temps1);
        mu_s=phi*sig_s.*temps2;
        S(Z(:,j),j)=randn(numz,1).*sqrt(sig_s)+mu_s;
        S(~Z(:,j),j)=randn(N-numz,1).*sqrt(1./alpha(j));
        X_k(:,Z(:,j)) = X_k(:,Z(:,j)) - sparse_mult(Index(:,Z(:,j)),A(:,j),S(Z(:,j),j));
    end
    info(it,X,A,S',Z',phi,alpha,'S');

    %sample alpha phi Pi
    ai=a0+sum(Z,1);
    bi=b0+N-sum(Z,1);
    Pi=betarnd(ai,bi);
    
    ci=c0+1/2*sum(sum(Index));
    di=d0+1/2*sum(sum((Index.*X-Index.*(A*(Z.*S)')).^2));
    phi=gamrnd(ci,1./di);
    
    ei=e0+K*N/2;
    fi=f0+1/2*sum(sum(S.^2));
    alpha(:)=gamrnd(ei,1./fi);
end
end

function DSsparse = sparse_mult(Yflag,D,S)
%DSsparse = (D*S').*Yflag;
%Use sparse multiplication to eliminate unnecessary computation
%Version 1: 12/01/2009
%Written by Mingyuan Zhou, Duke ECE, mz1@ee.duke.edu
[P,N]=size(Yflag);
DSsparse = sparse(1:P,1:P,double(D(:,1)))*Yflag*sparse(1:N,1:N,double(S(:,1)));
for k=2:size(D,2)
    DSsparse = DSsparse+ sparse(1:P,1:P,double(D(:,k)))*Yflag*sparse(1:N,1:N,double(S(:,k)));
end
end

function info(iter,Y,D,S,Z,ge,gs,msg)
rmse = sqrt(sum(sum((Y - D*(S.*Z)).^2)/numel(Y)));
psnr = -10*log10(rmse);
fprintf('%s %d\t%e\t%e\t%e\t%f\n', msg, iter, rmse, ge^-.5, mean(gs)^-.5, mean(sum(Z)));
end
