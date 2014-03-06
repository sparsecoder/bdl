classdef BPFA<handle
properties
    Y
    D,Z,S
    X,A,R,DTD
    pie, gs, ge
    P,N,K
    a,b,c,d,e,f
    
    sampleA=true, sampleD=true, sampleS=true, sampleZ=true
    verbose = true

end

methods
    function o = BPFA(Y, K)
        o.Y = Y;
        [o.P, o.N] = size(o.Y);
        o.K = K;
        o.a = o.K; o.b = o.N;
        o.c = 1e-6; o.d = 1e-6; o.e = 1e-6; o.f = 1e-6;
        
        o.D = o.Y(:,randperm(o.N,o.K));
        o.DTD = sum(o.D.^2);
        o.S = o.D'*o.Y;
        o.Z = o.S > mean(o.S(:)) + std(o.S(:));
        o.A = o.S.*o.Z;
        o.X = o.D*o.A;
        o.R = o.Y - o.X;
        o.pie = min(0.99999, sum(o.Z,2)/o.N);
        o.gs = 1;
        o.ge = 1e8;
    end
    
    function learn(o, T)
        if o.verbose
            fprintf('i (time): err \t\tgs\t\tge\t\tZiAvg\n');
            fprintf('%d: %e\t%e\t%e\t%e\n',0, o.err(), 1/o.gs, 1/o.ge, mean(sum(o.Z)/o.K) )
        end

        t=0;
        for i=1:T
            tic
            o.sample();
            t=t+toc;
            if o.verbose
                fprintf('%d (%f): %e\t%e\t%e\t%e\n',...
                  i, t, o.err(), 1/o.gs, 1/o.ge, mean(sum(o.Z)/o.K) )
            end
        end
    end
    
    function e = err(o)
        e = sqrt(sum(o.R(:).^2)/numel(o.X));
    end
    
    function sample(o)
        if o.sampleD, o.sample_D(); end

        if o.sampleA
            if o.sampleS, o.sample_S(); end
            if o.sampleZ, o.sample_Z(); end
        end
        
        o.sample_pie();
        o.sample_gs();
        o.sample_ge();
    end
    
    function sample_D(o)
        SigD = 1 ./( o.P + o.ge*sum(o.A.^2,2)*ones(1,o.P) );
        for k=1:o.K
            MuD = o.ge*SigD(k,:).*(o.A(k,:)*o.Xk(k)');
            d = mvnrnd(MuD', SigD(k,:)');

            o.X = o.X - (o.D(:,k) - d)*o.A(k,:);
            o.D(:,k) = d;
            o.DTD(k) = sum(d.^2);
            o.R = o.Y - o.X;
        end
    end
    
    function sample_S(o)
        for k = 1:o.K
            sigS = 1./(o.gs + o.ge*o.DTD(k));
            muS = o.ge*sigS*o.D(:,k)'*o.Xk(k);
            muS(o.Z(k,:)==0) = 0;
            
            s = randn(size(muS))*sigS + muS;

            o.X = o.X - o.D(:,k)*(o.A(k,:) - s.*o.Z(k,:));
            o.S(k,:) = s;
            o.R = o.Y - o.X;
        end
        o.A = o.S.*o.Z;
    end
    
    function sample_Z(o)
        for k = 1:o.K
            p1 = min(1e300, o.pie(k)*exp(...
                -o.ge/2*o.S(k,:).*(o.S(k,:)*o.DTD(k) - 2*o.D(:,k)'*o.Xk(k)) ));
            z = berrnd( p1./(p1 + 1 - o.pie(k)) );

            o.X = o.X - o.D(:,k)*(o.A(k,:) - o.S(k,:).*z);
            o.Z(k,:) = z;
            o.R = o.Y - o.X;
        end
        o.A = o.S.*o.Z;
    end
    
    function sample_pie(o)
        sumz = sum(o.Z,2);
        o.pie = min(0.99999, betarnd(o.a/o.K + sumz, o.b*(o.K-1)/o.K + o.N - sumz));
    end
    
    function sample_gs(o)
        o.gs = gamrnd(o.c + o.K*o.N/2, o.d + 2/sum(o.S(:).^2));
    end
    
    function sample_ge(o)
        o.ge = gamrnd(o.e + o.P*o.N, o.f + 2/norm(o.R, 'fro')^2);
    end
    
    function dtxk = DTXK(o,k,i)
        dtxk = o.D(:,k)'*o.Xk(k,i);
    end

    function xk = Xk(o,k)
        xk = o.R + o.D(:,k)*o.A(k,:);
    end
end

end

function r = berrnd(p)
    r = rand(size(p)) < p;
end

function r = mvnrnd(mu,sigma)
    r = randn(size(mu)).*sqrt(sigma) + mu; 
end

