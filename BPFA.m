classdef BPFA<handle
properties
    Y
    D,Z,S
    X,A,R,DTD
    pie, gs, ge
    P,N,K
    a,b,c,d,e,f
    
    sampleA, sampleD, sampleS, sampleZ

end

methods
    function o = BPFA(Y, K)
        o.Y = Y;
        [o.P, o.N] = size(o.Y);
        o.K = K;
        o.a = o.K; o.b = 1;
        o.c = 1e-6; o.d = 1e-6; o.e = 1e-6; o.f = 1e-6;
        
        o.D = normalize(o.Y(:,randperm(o.N,o.K)));
        o.S = o.D\o.Y;
        o.Z = o.S > mean(o.S(:)) - 1.8*std(o.S(:));
        o.Aup();
        o.pie = min(0.99999, sum(o.Z,2)/o.N);
        o.gs = 1e2;
        o.ge = 1e2;
        
        o.sampleA = true; o.sampleD = true;
        o.sampleS = true; o.sampleZ = true;
    end
    
    function learn(o, T)
        fprintf('i: err \t\tgs\t\tge\n');
        for i=1:T
            o.sample();
            fprintf('%d: %e\t%e\t%e\n',i, o.err(), 1/o.gs, 1/o.ge)
        end
    end
    
    function e = err(o)
        e = sqrt(sum(o.R(:).^2)/numel(o.X));
    end
    
    function sample(o)
        if o.sampleD
            o.sample_D();
            o.err()
        end
        
        if o.sampleA
            if o.sampleS, o.sample_S(); end
                o.err()
            if o.sampleZ, o.sample_Z(); end
                o.err()
        end
        
        o.sample_pie();
        o.sample_gs();
        o.sample_ge();
    end
    
    function Xup(o, k, a)
        if nargin < 2
            o.X = o.D*o.A;
            o.R = o.Y - o.X;
            o.DTD = sum(o.D.^2)';
        %else if nargin < 3
        %    o.X = o.D(:,k)*o.A;
        %    o.R = o.Y - o.X;
        %    o.DTD = sum(o.D.^2)';

        else
            o.X = o.X - o.D(:,k)*(o.A(k,:) - a);
            o.R = o.Y - o.X;
        end
    end
    function Aup(o,k,s,z)
        if nargin < 3
            o.A = o.S.*o.Z;
            o.Xup();
        else
            o.Xup(k,s.*z);
            o.S(k,:) = s;
            o.Z(k,:) = z;
            o.A(k,:) = s.*z; 
        end
    end
    
    function sample_D(o)
        SigD = 1 ./( o.P + o.ge*sum(o.A.^2,2)*ones(1,o.P) );
        for k=1:o.K
            MuD = o.ge*SigD(k,:).*(o.A(k,:)*o.Xk(k)');
            o.D(:,k) = mvnrnd(MuD', SigD(k,:))';
            o.Xup();
        end
    end
    
    function sample_S(o)
        %SigS = 1 ./( o.gs + o.ge*bsxfun(@times, o.Z.^2, o.DTD) );
        for k = 1:o.K
            sigS = 1./(o.gs + o.ge*o.DTD(k));
            dtxk = o.D(:,k)'*o.Xk(k);
            muS = o.ge*sigS*dtxk;
            muS(o.Z(k,:)==0) = 0;
            
            s = normrnd(muS', sigS)';

            o.Aup(k, s,o.Z(k,:));
        end
    end
    
    function sample_Z(o)
        for k = 1:o.K
            dtxk = o.D(:,k)'*o.Xk(k);
            for i=1:o.N
                p1(i) = min(1e300, o.pie(k)*exp(...
                    -o.ge/2*o.S(k,i)*(o.S(k,i)*o.DTD(k) - 2*dtxk(i)) ));
            end
            z = binornd( 1, p1./(p1 + 1 - o.pie(k)) );
            o.Aup(k, o.S(k,:),z);
        end
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

    function xk = Xk(o,k,i)
        if nargin < 3
            xk = o.R + o.D(:,k)*o.A(k,:);
        else
            xk = o.R(:,i) + o.D(:,k)*o.A(k,i);
        end
    end
end

end
