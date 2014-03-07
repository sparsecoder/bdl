classdef BPFA<handle
properties
    Y
    D,Z,S
    X,A,R,DTD
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
        o.a = o.K; o.b = o.N/8;
        o.c = 1e-6; o.d = 1e-6; o.e = 1e-6; o.f = 1e-6;
        
        o.D = o.Y(:,randperm(o.N,o.K));
        o.DTD = sum(o.D.^2);
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
        o.gs = 1e1;
        o.ge = 1e-2;
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
            o.print(); fprintf('\n');
        end
            fprintf('\n');

            o.sample_pie();
            o.sample_gs();
            o.sample_ge();
    end

    function sample_D(o,k)
       % SigD = 1 ./( o.P + o.ge*sum(o.A.^2,2)*ones(1,o.P) );
%        for k = o.ind
            SigD = 1 ./( o.P + o.ge*sum(o.A(k,:).^2)*ones(1,o.P) );
            xk = o.Y - o.D*(o.S.*o.Z) + o.D(:,k)*(o.S(k,:).*o.Z(k,:));
            MuD = o.ge*SigD.*(o.A(k,:)*xk');
%           % MuD = o.ge*SigD(k,:).*(o.A(k,:)*o.Xk(k)');
            d = mvnrnd(MuD', SigD');
%
%           % o.X = o.X - (o.D(:,k) - d)*o.A(k,:);
            o.D(:,k) = d;
            o.X = o.D*o.A;
            o.DTD(k) = sum(d.^2);
            o.R = o.Y - o.X;
%        end
    end
    
    function sample_S(o,k)
%        for k = o.ind
            sigS = 1./(o.gs + o.ge*o.Z(k,:)*sum(o.D(:,k).^2));
            muS = o.ge*sigS.*(o.D(:,k)'*o.Xk(k));
            muS(o.Z(k,:)==0) = 0;
            
            s = randn(size(muS)).*sigS + muS;

            o.X = o.X - o.D(:,k)*(o.A(k,:) - s.*o.Z(k,:));
            o.S(k,:) = s;
            o.R = o.Y - o.X;
%
%            for i = 1:o.N
%                if o.Z(k,i)
%                    xk = o.Y(:,i) - o.D*(o.S(:,i).*o.Z(:,i)) + o.D(:,k)*(o.S(k,i).*o.Z(k,i));
%                    dtxk = o.D(:,k)'*xk;
%                    sigS = 1./(o.gs + o.ge*sum(o.D(:,k).^2));
%                    muS = o.ge*sigS*dtxk;
%                else
%                    sigS = 1/o.gs;
%                    muS = 0;
%                end
%                s = randn()*sigS + muS;
%                
%                o.S(k,i) = s;
%                o.A(k,i) = o.Z(k,i)*s;
%                o.X = o.D*o.A;
%                o.R = o.Y - o.X;
%            end
%        end
%       % o.A = o.S.*o.Z;
    end
%    
    function sample_Z(o,k)
%        for k = o.ind
%            %dtxk = o.D(:,k)'*o.Xk(k);
%            for i = 1:o.N
%                xk = o.Y(:,i) - o.D*(o.S(:,i).*o.Z(:,i)) + o.D(:,k)*(o.S(k,i).*o.Z(k,i));
%                dtxk = o.D(:,k)'*xk;
%                dtd = sum(o.D(:,k).^2);
%                p1 = min(1e300, o.pie(k)*exp(...
%                    -o.ge/2*o.S(k,i).*(o.S(k,i)*dtd - 2*dtxk) ));
%                z = berrnd( p1./(1 - o.pie(k) + p1) );
%
%                %o.X(:,i) = o.X(:,i) - o.D(:,k)*(o.A(k,i) - o.S(k,i).*z);
%                %o.Z(k,i) = z;
%                %o.R(:,i) = o.Y(:,i) - o.X(:,i);
%                
%                o.Z(k,i) = z;
%                o.A(k,i) = o.S(k,i)*z;
%                o.X = o.D*o.A;
%                o.R = o.Y - o.X;
%            end
%        end
%       % o.A = o.S.*o.Z;
        dtxk = o.D(:,k)'*o.Xk(k);
        p1 = min(1e300, o.pie(k)*exp(...
             -o.ge/2*o.S(k,:).*(o.S(k,:)*sum(o.D(:,k).^2) - 2*dtxk) ));
        z = berrnd( p1./(1 - o.pie(k) + p1) );
        o.X = o.X - o.D(:,k)*(o.A(k,:) - o.S(k,:).*z);
        o.Z(k,:) = z;
        o.R = o.Y - o.X;

    end
    
    function xk = Xk(o,k)
        xk = o.R + o.D(:,k)*o.A(k,:);
    end

    function sample_pie(o)
        sumz = sum(o.Z,2);
        o.pie = min(0.99, betarnd(o.a/o.K + sumz, o.b*(o.K-1)/o.K + o.N - sumz));
    end
    
    function sample_gs(o)
        o.gs = max(1e-5, gamrnd(o.c + o.K*o.N/2, o.N/(o.d + 0.5*sum(o.S(:).^2)) ));
    end
    
    function sample_ge(o)
        o.ge = max(1e-4, gamrnd(o.e + o.P*o.N/2, o.N/(o.f + 0.5*sum((o.Y(:)-o.X(:)).^2)) ));
    end
    
    function learn(o, T)
        if o.verbose
            s = ' ';
            fprintf('iter:\ttime\tpsnr%*cpsnr0%*cgs^-.5%*cge^-.5%*cZfill\n',...
                8,s,6,s,6,s,6,s);
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
        [o.psnr(o.R), o.psnr(o.X0-o.X), o.gs^-.5, o.ge^-.5, nnz(o.Z)/numel(o.Z)] );
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

