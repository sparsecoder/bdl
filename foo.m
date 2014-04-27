function [D, S, Z, ge, gs, Pi] = foo(Y,K, iter, D, ge, gs, Pi, S, Z)
[P,N] = size(Y);

sampleD = true; defD = randn(P,K);
sampleS = true; defS = randn(K,N);
sampleZ = true; defZ = false(K,N);
sampleGE = true; defGE = 1;
sampleGS = true; defGS = ones(K,1);
samplePi = true; defPi = ones(K,1);

    function [sample, arg] = checkArg(arg, default)
        sample = false;
        if isempty(arg)
            arg = default;
            sample = true;
        end
    end

if nargin>=4, [sampleD, D] = checkArg(D, defD); K=size(D,2); else D = defD;  end
if nargin>=5, [sampleGE, ge] = checkArg(ge, defGE); else ge = defGE; end
if nargin>=6, [sampleGS, gs] = checkArg(gs, defGS); else gs = defGS; end
if nargin>=7, [samplePi, Pi] = checkArg(Pi, defPi); else Pi = defPi; end
if nargin>=8, [sampleS, S] = checkArg(S, defS); else S = defS; end
if nargin>=9, [sampleZ, Z] = checkArg(Z, defZ); else Z = defZ; end


gea0 = 1e-6;
geb0 = 1e-6;
gsa0 = 1e-6;
gsb0 = 1e-6;
pia0 = 1;
pib0 = N/8;

R = Y - D*(S.*Z);
for iter = 1:iter
    if sampleD
        for k = 1:K
            R(:,Z(k,:)) = R(:,Z(k,:)) + D(:,k)*S(k,Z(k,:));        
            sig = 1./(P + ge*sum(S(k,Z(k,:))'.^2));
            mu = ge*sig.*(R(:,Z(k,:))*S(k,Z(k,:))');
            D(:,k) = mu + sqrt(sig).*randn(P,1);
            R(:,Z(k,:)) = R(:,Z(k,:)) - D(:,k)*S(k,Z(k,:));
        end
        info(iter,R,Z,ge,gs,Pi,'D');
    end
    
    if sampleZ
        for k = 1:K
            R(:,Z(k,:)) = R(:,Z(k,:)) + D(:,k)*S(k,Z(k,:));
            dtd = sum(D(:,k).^2);
            dtxk = D(:,k)'*R;
            p1 = Pi(k)*exp( -0.5*ge*(S(k,:).^2*dtd - 2*S(k,:).*dtxk) );
            Z(k,:) = rand(1,N) > ( (1 - Pi(k))./(1 - Pi(k) + p1) );
            R(:,Z(k,:)) = R(:,Z(k,:)) - D(:,k)*S(k,Z(k,:));
        end
        info(iter,R,Z,ge,gs,Pi,'Z');
    end
    
    if sampleS
        for k = 1:K
            R(:,Z(k,:)) = R(:,Z(k,:)) + D(:,k)*S(k,Z(k,:));
            nz = nnz(Z(k,:));
            dtd = sum(D(:,k).^2);
            dtxk = D(:,k)'*R(:,Z(k,:));
            sig = 1./(gs(k)*ones(1,nz) + ge*dtd);
            mu = ge*sig.*dtxk;
            S(k,Z(k,:)) = mu + sqrt(sig).*randn(1,nz);
            S(k,~Z(k,:)) = sqrt(1/gs(k)).*randn(1,N-nz);
            R(:,Z(k,:)) = R(:,Z(k,:)) - D(:,k)*S(k,Z(k,:));
        end
        info(iter,R,Z,ge,gs,Pi,'S');
    end
    
    R = Y - D*(S.*Z);
    
    if sampleGE
        geai = gea0 + 0.5*numel(Y);
        gebi = geb0 + 0.5*sum(R(:).^2);
        ge = gamrnd(geai,1/gebi);
    end
    
    if sampleGS
        gsai = gsa0 + 0.5*numel(S);
        gsbi = gsb0 + 0.5*sum(sum(S.^2,2));% + 0.5*(N - sum(Z,2))./gs;
        gs(:) = gamrnd(gsai, 1./gsbi);
    end
    
    if samplePi
        piai = pia0 + sum(Z,2);
        pibi = pib0 + N - sum(Z,2);
        Pi = betarnd(piai,pibi);
    end
    info(iter,R,Z,ge,gs,Pi,'P');
end

end

function info(iter,R,Z,ge,gs,Pi,msg)
rmse = sqrt(sum(sum(R.^2)/numel(R)));
psnr = -10*log10(rmse);
fprintf('%s %d\t%e\t%e\t%e\t%f\t%f\n', msg, iter, rmse, ge^-.5, mean(gs)^-.5, mean(sum(Z)), mean(Pi));
end
