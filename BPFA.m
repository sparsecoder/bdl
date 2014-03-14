classdef BPFA<handle
properties
    Y
    D,Z,S
    X,A,R
    X0
    pie, gs, ge
    P,N,K
    a,b,c,d,e,f
    ind
    
    sampleA=true, sampleD=true, sampleS=true, sampleZ=true
    verbose = true

end

methods
    function o = BPFA(Y, X0, K)
        o.Y = Y;
        o.X0 = X0;
        [o.P, o.N] = size(o.Y);
        o.K = K;
        o.a = o.K; o.b = 1;
        o.c = 1e-6; o.d = 1e-6; o.e = 1e-6; o.f = 1e-6;
        
        o.D = o.Y(:,randperm(o.N,2*o.K));
        [~,ind] = sort(std(o.D),'descend');
        o.D = o.D(:,ind(1:o.K));
        displayPatches(o.D); drawnow;

        o.S = o.D'*o.Y;
        o.Z = o.S > mean(o.S(:)) + 2.5*std(o.S(:));
        o.init(o)
    end

    function init(o, mats)
        if isfield(mats, 'D'), o.D = mats.D; end
        if isfield(mats, 'S'), o.S = mats.S; end
        if isfield(mats, 'Z'), o.Z = mats.Z; end
        o.Z = sparse(logical(o.Z));
          
        o.pie = min(0.99999, sum(o.Z,2)/o.N);
        o.gs = 1/cov(o.Z(:).*o.S(:));
        Err2 = (o.Y - o.D*(o.S.*o.Z)).^2;
        o.ge = 1;
        o.A = o.S.*o.Z;
        o.X = o.D*(o.A);
        o.R = o.Y - o.X;
        o.check();
    end
    
    function check(o)
        if any(isnan([o.D(:); o.S(:)])) %; o.Z(:)])) z is boolean...
            error('nan found');
        end
    end

    function sample(o)
        k = randperm(o.K);
        if o.sampleD, o.sample_D(k); end

        if o.sampleA
            if o.sampleS, o.sample_S(k); end
            if o.sampleZ, o.sample_Z(k); end
        end

        o.sample_pie();
        o.sample_gs();
        o.sample_ge();
    end

    function sample_D(o,K)
        %sig = bsxfun(@times, o.P + o.ge*sum(o.A'.^2), ones(o.P,o.K)).^-1;
        %sig = (o.P + o.ge*sum(o.A.^2,2)').^-1;
        for k=K
            if nnz(o.Z(k,:))
                xk = o.Y(:,o.Z(k,:)) - o.X(:,o.Z(k,:)) + o.D(:,k)*o.S(k,o.Z(k,:));
                sig = 1./(o.P + o.ge*sum(o.S(k,o.Z(k,:)).^2));
                mu = o.ge*sig.*(xk*o.S(k,o.Z(k,:))');
            else
                sig = 1/o.P;
                mu = zeros(o.P,1);
            end
            d = mvnrnd(mu,sig);

            if nnz(o.Z(k,:))
                o.X(:,o.Z(k,:)) = o.X(:,o.Z(k,:)) + (d - o.D(:,k))*o.S(k,o.Z(k,:));
            end
            o.D(:,k) = d;
            %if nnz(o.Z(k,:)), o.R(:,o.Z(k,:)) = o.R(:,o.Z(k,:)) - o.D(:,k)*o.S(k,o.Z(k,:)); end
        end
    end
    
    function sample_S(o,K)
        dtd = sum(o.D.^2)';
        sig = bsxfun(@times, o.gs + o.ge*o.Z, dtd).^-1;
        for k=K
            xk = o.Y(:,o.Z(k,:)) - o.X(:,o.Z(k,:)) + o.D(:,k)*o.S(k,o.Z(k,:));
            dtxk = o.D(:,k)'*xk;
            mu = zeros(1,o.N);
            mu(o.Z(k,:)) = o.ge*sig(k,o.Z(k,:)).*dtxk;
            s = randn(1,o.N).*sqrt(sig(k,:)) + mu;
            if nnz(o.Z(k,:))
                o.X(:,o.Z(k,:)) = o.X(:,o.Z(k,:))...
                     + o.D(:,k)*(s(o.Z(k,:)) - o.S(k,o.Z(k,:)));
            end
            o.S(k,:) = s;
        end
    end
    
    function sample_Z(o,K)
        dtd = sum(o.D.^2)';
        sdtd = bsxfun(@times, o.S.^2, dtd);
        for k=K
            xk = o.Y - o.X + o.D(:,k)*(o.S(k,:).*o.Z(k,:));
            dtxk = o.D(:,k)'*xk;
            p1 = o.pie(k)*exp(-0.5*o.ge*(sdtd(k,:) - 2*o.S(k,:).*dtxk));
            z = berrnd(p1./(1 - o.pie(k) + p1));

            o.X = o.X + o.D(:,k)*(o.S(k,:).*z - o.S(k,:).*o.Z(k,:));
            o.Z(k,:) = z;
        end
    end
    
    function sample_ge(o)
        Err2 = (o.Y - o.D*(o.S.*o.Z)).^2;
        o.ge = gamrnd(o.e + o.P*o.N/2, 1/(o.f + 0.5*sum(Err2(:))) );
    end
    
    function sample_gs(o)
        %o.gs = gamrnd(o.c + o.K*o.N, 1/(o.d + 0.5*sum(o.S(:).^2)) );
        %o.gs = gamrnd(o.c + nnz(o.Z), 1/(o.d + 0.5*sum((o.Z(:).*o.S(:)).^2)) );
        o.gs = 1;
    end
    
    function sample_pie(o)
        sumz = sum(o.Z,2);
        o.pie = betarnd(o.a/o.K + sumz, o.b*(o.K-1)/o.K + o.N - sumz);
    end
    
    function learn(o, T)
        if o.verbose
            s = ' ';
            fprintf('iter:\ttime\tpsnr%*cpsnr0%*csds%*csde%*cZfill\n',...
                8,s,6,s,8,s,8,s);
            fprintf('%4d: %6.0fs\t', 0, 0);
            o.print();
        end

        t=0;
        for i=1:T
            tic
            o.sample();
            t=t+toc;
            if o.verbose
                fprintf('%4d: %6.0fs\t',i, t);
                o.print();
                [~,ind] = sort(sum(o.Z,2),'descend');
                displayPatches(o.D(:,ind)); drawnow;
            end
        end
    end
    
    function print(o)
        fprintf('%3.4e ',...
        [o.rms(o.Y-o.D*(o.S.*o.Z)), o.rms(o.X0-o.D*(o.S.*o.Z)),...
          o.gs^-0.5, o.ge^-0.5, nnz(o.Z)/numel(o.Z)] );
        fprintf('\n');
    end

    function e = psnr(o,x)
        e = 20*log10(1/o.rms(x));
    end
    function e = rms(o, x)
        e = sqrt(sum(sum(x.^2))/numel(x)); 
    end
    
end

end

function r = berrnd(p)
    r = rand(size(p)) < p;
end

function r = mvnrnd(mu,sigma)
    r = randn(size(mu)).*sqrt(sigma) + mu; 
end

