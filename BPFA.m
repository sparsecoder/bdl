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
        
        o.D = o.Y(:,randperm(o.N,o.K));
        o.S = o.D'*o.Y;
        o.Z = o.S > mean(o.S(:)) + std(o.S(:));
        o.init(o)
    end

    function init(o, mats)
        if isfield(mats, 'D'), o.D = mats.D; end
        if isfield(mats, 'S'), o.S = mats.S; end
        if isfield(mats, 'Z'), o.Z = mats.Z; end
             
        o.A = o.S.*o.Z;
        o.X = o.D*o.A;
        o.R = o.Y - o.X;
        o.pie = min(0.99999, sum(o.Z,2)/o.N);
        o.gs = 1e-3;
        o.ge = 0.1;
    end
    
    function sample(o)
        for k=1:o.K; % randperm(o.K);
            if o.sampleD, o.sample_D(k); end
        end
        o.print(); fprintf('\n');

        if o.sampleA
            for k=1:o.K; % randperm(o.K);
                if o.sampleS, o.sample_S(k); end
            end
            o.print(); fprintf('\n');
            for k=1:o.K; % randperm(o.K);
                if o.sampleZ, o.sample_Z(k); end
            end
            %o.print(); fprintf('\n');
        end
        %end
        fprintf('\n');

        o.sample_pie();
        o.sample_gs();
        o.sample_ge();
    end

    function check(o)
        if any(isnan([o.D(:); o.S(:); o.Z(:)]))
            error('nan found');
        end
    end

    function sample_D(o,k)
        o.check();
        A = o.S.*o.Z;
        sig = (o.P + o.ge*sum(A(k,:).^2)*ones(o.P,1)).^-1;
        xk = o.Y - o.D*A + o.D(:,k)*A(k,:);
        mu = o.ge*sig.*(xk*A(k,:)');
        o.D(:,k) = mvnrnd(mu,sig);

        o.X = o.D*A;
        o.R = o.Y - o.X;
    end
    
    function sample_S(o,k)
        o.check();
        dtd = sum(o.D(:,k).^2);
        sig = (o.gs + o.ge*o.Z(k,:)*dtd).^-1;
        mu = zeros(1,o.N);
        I = find(o.Z(k,:));
        A = o.S(:,I).*o.Z(:,I);
        xk = o.Y(:,I) - o.D*A + o.D(:,k)*A(k,:);
        dtxk = o.D(:,k)'*xk;
        mu(I) = o.ge*sig(I).*dtxk;
        o.S(k,:) = randn(1,o.N).*sig + mu;

        o.X = o.D*(o.S.*o.Z);
        o.R = o.Y - o.X;
    end
    
    function sample_Z(o,k)
        o.check();
        dtd = sum(o.D(:,k).^2);
        A = o.S.*o.Z;
        xk = o.Y - o.D*A + o.D(:,k)*A(k,:);
        dtxk = o.D(:,k)'*xk;
        p1 = o.pie(k)*exp(-o.ge/2*o.S(k,:).*(o.S(k,:)*dtd - 2*dtxk));
        o.Z(k,:) = berrnd(p1./(1 - o.pie(k) + p1));

        o.X = o.D*(o.S.*o.Z);
        o.R = o.Y - o.X;
    end
    
    function sample_ge(o)
        o.ge = gamrnd(o.e + o.P*o.N/2, 1/(o.f + 0.5*sum(o.R(:).^2)));
    end
    
    function sample_gs(o)
        o.gs = gamrnd(o.c + o.K*o.N/2, 1/(o.d + 0.5*sum(o.S(:).^2)));
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
            end
        end
    end
    
    function print(o)
        fprintf('%3.4e ',...
        [o.psnr(o.R), o.psnr(o.X0-o.X), o.gs^-0.5, o.ge^-0.5, nnz(o.Z)/numel(o.Z)] );
        fprintf('\n');
    end

    function e = psnr(o,x)
        e = 20*log10(255/o.rms(x));
    end
    function e = rms(o, x)
        e = sqrt(mean(sum(x(:).^2))); 
    end
    
end

end

function r = berrnd(p)
    r = rand(size(p)) < p;
end

function r = mvnrnd(mu,sigma)
    r = randn(size(mu)).*sqrt(sigma) + mu; 
end

