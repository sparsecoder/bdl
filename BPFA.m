classdef BPFA<handle
properties
    Y
    D,Z,S
    X,A,R
    pie, gs, ge
    P,N,K
    a,b,c,d,e,f
end

methods
    function o = BPFA(Y)
        o.Y = Y;
        [o.P, o.N] = size(o.Y);
        o.K = 16;
        o.a = 1e-6; o.b = 1e-6; o.c = 1e-6; o.d = 1e-6; o.e = 1e-6; o.f = 1e-6;
        
        o.D = o.Y(:,randperm(o.N,o.K));
        o.S = o.D'*o.Y;
        o.Z = o.S-mean(o.S(:)) > 0;
        o.A = o.S.*o.Z;
        o.X = o.D*o.A;
        o.R = o.Y - o.X;
        o.pie = ones(o.K,1)/o.K;
        o.gs = 1;
        o.ge = 1;
        
        for i=1:100
            o.sample();
            fprintf('%d: %f\n',i, norm(o.R, 'fro'))
        end
        
    end
    
    function sample(o)
        o.sample_D();
        o.sample_Z();
        o.sample_S();
        o.sample_pie();
        o.sample_gs();
        o.sample_ge();
        
        o.A = o.S.*o.Z;
        o.X = o.D*o.A;
        o.R = o.Y - o.X;
    end
    
    function sample_D(o)
        SigmaD = zeros(1,o.P,o.K) + o.P;
        MuD = zeros(size(o.D));
        for k=1:o.K
            SigmaD(:,:,k) = SigmaD(:,:,k)...
                + o.ge*sum(o.A(:,k).^2)*ones(1,o.P);
            
            MuD(:,k) = o.ge*SigmaD(:,:,k).*(o.A(k,:)*o.Xk(k)');
        end
        
        SigmaD = 1./SigmaD;
        o.D = mvnrnd(MuD',SigmaD)';
    end
    
    function sample_Z(o)
        dtd = o.DTD();
        p1 = zeros(size(o.Z));
        for k=1:o.K
            p1(k,:) = o.S(k,:)*dtd(k) - 2*o.DTXK(k);
        end
        p1 = bsxfun(@times, exp( -o.ge/2*(o.S.*p1) ), o.pie);
        den = bsxfun(@plus, p1, 1-o.pie);
        
        o.Z = binornd(1,p1./den);
    end
    
    function sample_S(o)
        SigmaS = 1 ./( o.gs + o.ge* bsxfun(@times, o.Z.^2, o.DTD()) );
        MuS = zeros(size(SigmaS));
        for k=1:o.K
            MuS(k,:) = o.ge*SigmaS(k,:).*o.Z(k,:).*o.DTXK(k);
        end

        o.S = normrnd(MuS,SigmaS);
    end
    
    function sample_pie(o)
        sumz = sum(o.Z);
        o.pie = betarnd(o.a/o.K + sumz, o.b*(o.K-1)/o.K + o.N - sumz);
    end
    
    function sample_gs(o)
        o.gs = gamrnd(o.c + o.K*o.N/2, o.d + norm(o.S, 'fro')^2/2);
    end
    
    function sample_ge(o)
        o.ge = gamrnd(o.e, o.f + norm(o.R, 'fro')^2/2);
    end
    
    function dtd = DTD(o)
        dtd = dot(o.D,o.D)';
    end
    
    function dtxk = DTXK(o,k)
        dtxk = o.D(:,k)'*o.Xk(k);
    end
    
    function sts = STS(o)
        sts = dot(o.S,o.S);
    end

    function xk = Xk(o,k)
        xk = o.R + o.D(:,k)*o.A(k,:); %!!! this is 1x1
    end
end

end