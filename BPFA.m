classdef BPFA < handle 

properties
    Y,X0
    X,R
    D,S,Z
    ge,gs,Pi
    
    K

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
        if isfield(args, 'X'), o.X0 = args.X; else, o.X0 = args.Y; end
        if isfield(args, 'D'), o.D = args.D; end
        if isfield(args, 'S'), o.S = args.S; end
        if isfield(args, 'Z'), o.Z = args.Z; end
        o.Z = logical(o.Z);

        if isfield(args, 'ge'), o.ge = args.ge; else, o.ge = 1; end
        if isfield(args, 'gs'), o.gs = args.gs; else o.gs = 1; end
        o.pia = o.K;
        o.pib = 1;

        if isfield(args, 'inferA'), o.inferA = args.inferA; end
        if isfield(args, 'inferD'), o.inferD = args.inferD; end
        if isfield(args, 'inferS'), o.inferS = args.inferS; end
        if isfield(args, 'inferZ'), o.inferZ = args.inferZ; end
        if isfield(args, 'inferGe'), o.inferGe = args.inferGe; end
        if isfield(args, 'inferGs'), o.inferGs = args.inferGs; end
    end

    function learn(o)
        for i = 1:20
            o.R = o.Y - o.D*(o.S.*o.Z);
            if o.inferGe, o.ge = o.sample_ge(); end
            if o.inferGs, o.gs = o.sample_gs(); end
            o.Pi = o.sample_pi();
        end


    end

    function ge = sample_ge(o)
        ge = gamrnd(o.gea + numel(o.Y)/2, 1/(o.geb + 0.5*sum(o.R(:).^2)) );
    end

    function gs = sample_gs(o)
        gs = gamrnd(o.gsa + numel(o.S), 1/(o.gsb + 0.5*(sum(o.S(:).^2) + numel(o.S) - nnz(o.Z)/o.gs)) );
    end

    function Pi = sample_pi(o)
        sumz = sum(o.Z,2);
        Pi = betarnd(o.pia/o.K + sumz, o.pib*(o.K-1)/o.K + size(o.Y,2) - sumz);
    end

end

end

