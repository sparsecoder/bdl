classdef BPFA<handle
properties
    Y
    D,Z,S
    X,A,R,DTD
    pie, gs, ge
    P,N,K
    a,b,c,d,e,f
    
    sampleA, sampleD
end

methods
    function o = BPFA(Y, K)
        o.Y = Y;
        [o.P, o.N] = size(o.Y);
        o.K = K;
        o.a = o.N*o.K; o.b = o.a*(1/0.7 + 0.70);
        o.c = 1e-6; o.d = 1e-6; o.e = 1e-6; o.f = 1e-6;
        
        o.D = normalize(o.Y(:,randperm(o.N,o.K)));
        o.S = o.D\o.Y;
        o.Z = o.S > mean(o.S(:)) - 1.8*std(o.S(:));
        o.Aup();
        o.pie = ones(o.K,1)/o.K;
        o.gs = 10^3;
        o.ge = 10^3;
        
        o.sampleA = true; o.sampleD = true;
    end
    
    function learn(o, T)
        fprintf('i: err \t\tgs\t\tge\n');
        for i=1:T
            o.sample();
            fprintf('%d: %f\t%f\t%f\n',i, o.err(), 1/o.gs, 1/o.ge)
        end
    end
    
    function e = err(o)
        e = sqrt(sum(o.R(:).^2)/numel(o.X));
    end
    
    function sample(o)
        if o.sampleD
            o.sample_D();
        end
        
        if o.sampleA
            o.sample_Z();
            o.Aup();
            o.sample_S();
        end
        o.Aup();
        
        o.sample_pie();
        o.sample_gs();
        o.sample_ge();
    end
    
    function Xup(o)
        o.X = o.D*o.A;
        o.R = o.Y - o.X;
        o.DTD = dot(o.D,o.D)';
    end
    function Aup(o)
        o.A = o.S.*o.Z;
        o.Xup();
    end
    
    function sample_D(o)
        SigD = 1 ./( o.P + o.ge*sum(o.A.^2,2)*ones(1,o.P) );
        MuD = zeros(size(o.D));
        for k=1:o.K
            MuD(:,k) = o.ge*SigD(k,:).*(o.A(k,:)*o.Xk(k)');
        end
        
        o.D = mvnrnd(MuD', reshape(SigD, [1 size(SigD')]))';
    end
    
    function sample_Z(o)
%         p1_1 = zeros(size(o.Z));
%         for k=1:o.K
%             p1_1(k,:) = o.S(k,:).*(o.S(k,:)*dtd(k) - 2*o.DTXK(k));
%         end
%         p1_1 = bsxfun(@times, exp(-o.ge/2*p1_1), o.pie);  %/median(p1(:))
        
        
        p1 = zeros(size(o.Z));
        for k = 1:o.K
            dtxk = o.DTXK(k);
            for i=1:o.N
                p1(k,i) = o.pie(k)*exp(...
                    -o.ge/2*(o.S(k,i).^2*o.DTD(k) - 2*o.S(k,i)*dtxk(i)) );
            end
        end
        
        p1(p1==inf) = 1;
        den = bsxfun(@plus, p1, 1-o.pie);
        
        o.Z = binornd(1,p1./den);
    end
    
    function sample_S(o)
        SigS = 1 ./( o.gs + o.ge*bsxfun(@times, o.Z.^2, o.DTD) );
%         MuS1 = zeros(size(o.S));
%         for k=1:o.K
%             MuS1(k,:) = o.ge*SigS1(k,:).*o.Z(k,:).*o.DTXK(k);
%         end

        MuS = zeros(size(o.S));
        for k = 1:o.K
            dtxk = o.DTXK(k);
            for i=1:o.N
                MuS(k,i) = o.ge*SigS(k,i)*dtxk(i);
            end
        end

        o.S = normrnd(MuS,SigS);
    end
    
    function sample_pie(o)
        sumz = sum(o.Z);
        o.pie = betarnd(o.a/o.K + sumz, o.b*(o.K-1)/o.K + o.N - sumz);
    end
    
    function sample_gs(o)
        o.gs = gamrnd(o.c + o.K*o.N/2, 2/(o.d + sum(dot(o.S,o.S))));
    end
    
    function sample_ge(o)
        o.ge = gamrnd(o.e + o.P*o.N, 2/(o.f + norm(o.R, 'fro')^2));
    end
    
    function dtxk = DTXK(o,k)
        dtxk = o.D(:,k)'*o.Xk(k);
    end

    function xk = Xk(o,k)
        xk = o.R + o.D(:,k)*o.A(k,:);
    end
end

end