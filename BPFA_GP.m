function [D, S, Z, phi, alpha] = BPFA_GP(X,K,iter, sampleD, D, phi, alpha)
[P N] = size(X);

c0 = 1e-6;
d0 = 1e-6;
e0 = 1e-6;
f0 = 1e-6;
a0 = 1;
b0 = N/8;

if nargin<4, sampleD = true; end
if sampleD, D = randn(P,K); else K=size(D,2); end
if nargin<6, phi = 1; samplePhi = true; else samplePhi = false; end
if nargin<7, alpha = ones(1,K); sampleAlpha = true; else sampleAlpha = false; end
S = randn(N,K);
Z = false(N,K);
Pi = ones(K,1);

R = X - D*(Z.*S)';
info(0,R,Z',phi,alpha,'0');
for it = 1:iter
    if sampleD
        for k = 1:K
            R(:,Z(:,k)) = R(:,Z(:,k)) + D(:,k)*S(Z(:,k),k)';
            sig = 1./(P + phi*sum(S(Z(:,k),k).^2));
            mu = phi*sig.*(R(:,Z(:,k))*S(Z(:,k),k));
            D(:,k) = mu + sqrt(sig).*randn(P,1);
            R(:,Z(:,k)) = R(:,Z(:,k)) - D(:,k)*S(Z(:,k),k)';
        end
        info(it,R,Z',phi,alpha,'D');
    end
    
    for k = 1:K
        R(:,Z(:,k)) = R(:,Z(:,k)) + D(:,k)*S(Z(:,k),k)';
        dtd = sum(D(:,k).^2);
        dtxk = R'*D(:,k);
        p1 = Pi(k)*exp( -0.5*phi*(S(:,k).^2.*dtd - 2*S(:,k).*dtxk) );
        Z(:,k) = rand(N,1)>((1-Pi(k))./(1-Pi(k)+p1));
        R(:,Z(:,k)) = R(:,Z(:,k)) - D(:,k)*S(Z(:,k),k)';
    end
    info(it,R,Z',phi,alpha,'Z');
    
    for k = 1:K
        R(:,Z(:,k)) = R(:,Z(:,k)) + D(:,k)*S(Z(:,k),k)';
        numz = nnz(Z(:,k));
        dtd = sum(D(:,k).^2);
        dtxk = R(:,Z(:,k))'*D(:,k);
        sig = 1./(alpha(k)*ones(numz,1) + phi*dtd);
        mu = phi*sig.*dtxk;
        S(Z(:,k),k) = mu + sqrt(sig).*randn(numz,1);
        S(~Z(:,k),k) = sqrt(1/alpha(k)).*randn(N-numz,1);
        R(:,Z(:,k)) = R(:,Z(:,k)) - D(:,k)*S(Z(:,k),k)';
    end
    info(it,R,Z',phi,alpha,'S');
    
    R = X - D*(Z.*S)';
    
    if samplePhi
        ci = c0+1/2*numel(X);
        di = d0+1/2*sum(sum(R.^2));
        phi = gamrnd(ci,1./di);
    end
    
    if sampleAlpha
        ei = e0+numel(Z)/2;
        fi = f0+1/2*sum(sum(S.^2));
        alpha(:) = gamrnd(ei,1./fi);
    end
    
    ai = a0+sum(Z,1);
    bi = b0+N-sum(Z,1);
    Pi = betarnd(ai,bi);
end
end

function info(iter,R,Z,ge,gs,Pi,msg)
rmse = sqrt(sum(sum(R.^2)/numel(R)));
psnr = -10*log10(rmse);
fprintf('%s %d\t%e\t%e\t%e\t%f\t%f\n', msg, iter, rmse, ge^-.5, mean(gs)^-.5, mean(sum(Z)), mean(Pi));
end

