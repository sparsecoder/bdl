classdef BPFA < handle 

properties
    Y,X0
    X,R
    D,S,Z
    ge,gs,Pi
    
    K,P,N

    gea = 1e-6, geb = 1e-6
    gsa = 1e-6, gsb = 1e-6
    pia,pib
    
    inferA = true
    inferD = true
    inferS = true, inferZ = true
    inferGe = true, inferGs = true
end

methods
    function o = BPFA(args)
        o.Y = args.Y;
        o.K = args.K;
        [o.P,o.N] = size(o.Y);
        if isfield(args, 'X'), o.X0 = args.X; else, o.X0 = args.Y; end
        if isfield(args, 'D'), o.D = args.D; else, o.D = o.Y(:,randperm(o.N,o.K)); end
        if isfield(args, 'S'), o.S = args.S; end
        if isfield(args, 'Z'), o.Z = args.Z; end
        o.Z = logical(o.Z);
        o.X = o.D*(o.S.*o.Z);

        if isfield(args, 'ge'), o.ge = args.ge; else, o.ge = 1; end
        if isfield(args, 'gs'), o.gs = args.gs; else o.gs = 1; end
        o.pia = o.K;
        o.pib = 1;

        if isfield(args, 'inferA'), o.inferA = args.inferA; end
        if isfield(args, 'inferD'), o.inferD = args.inferD; end
        if isfield(args, 'inferS'), o.inferS = args.inferS; end
        if isfield(args, 'inferZ'), o.inferZ = args.inferZ; end
        o.inferS = o.inferA; o.inferZ = o.inferA;
        if isfield(args, 'inferGe'), o.inferGe = args.inferGe; end
        if isfield(args, 'inferGs'), o.inferGs = args.inferGs; end
    end

    function learn(o)
        for i = 1:20
            if o.inferD, o.sample_D(); end
            o.R = o.Y - o.X;
            if o.inferGe, o.ge = o.sample_ge(); end
            if o.inferGs, o.gs = o.sample_gs(); end
            o.Pi = o.sample_pi();

            %o.info();
        end


    end

    function sample_D(o)
        for k=1:o.K
            if nnz(o.Z(k,:))
                xk = o.Y(:,o.Z(k,:)) - o.X(:,o.Z(k,:)) + o.D(:,k)*o.S(k,o.Z(k,:));
                sig = 1./(o.P + o.ge*sum(o.S(k,o.Z(k,:)).^2));
                mu = o.ge*sig.*(xk*o.S(k,o.Z(k,:))');
            else
                sig = 1/o.P;
                mu = zeros(o.P,1);
            end
            d = mu + sqrt(sig).*randn(size(mu));

            if nnz(o.Z(k,:))
                o.X(:,o.Z(k,:)) = o.X(:,o.Z(k,:)) + (d - o.D(:,k))*o.S(k,o.Z(k,:));
            end
            o.D(:,k) = d;


%            xk = o.Y - o.D*(o.S.*o.Z) + o.D(:,k)*(o.S(k,:).*o.Z(k,:));
%            sig = 1./(o.P + o.ge*sum((o.S(k,:).*o.Z(k,:)).^2));
%            mu = o.ge*sig.*(xk*(o.S(k,:).*o.Z(k,:))');
%            
%            d = mu + sqrt(sig)*randn(size(mu));
%
%            o.D(:,k) = d;
%            o.X = o.D*(o.S.*o.Z); %+ o.X + (d - o.D(:,k))*(o.S(k,:).*o.Z(k,:));
            
        end
    end

    function ge = sample_ge(o)
        ge = gamrnd(o.gea + numel(o.Y)/2, 1/(o.geb + 0.5*sum(o.R(:).^2)) );
    end

    function gs = sample_gs(o)
        gs = gamrnd(o.gsa + numel(o.S)/2, 1/(o.gsb + 0.5*(sum(o.S(:).^2))));% + numel(o.S) - nnz(o.Z)/o.gs)) );
    end

    function Pi = sample_pi(o)
        sumz = sum(o.Z,2);
        Pi = betarnd(o.pia/o.K + sumz, o.pib*(o.K-1)/o.K + o.N - sumz);
    end

    function info(o)
        R0 = o.X0 - o.X;
        fprintf('%f\t%f\t%f\n',sqrt(mean(R0(:).^2)), sqrt(mean(o.R(:).^2)), o.ge);
        
    end

end

end

