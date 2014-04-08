classdef BPFATest < matlab.unittest.TestCase

properties
    Y,X
    D,S,Z
    ge,gs

    K

    inferA = false
    inferD = false, inferS = false, inferZ = false
    inferGe = false,inferGs = false

    verbose = true
end

methods (TestMethodSetup)
    function setup(o)
        warning('off', 'MATLAB:structOnObject');
        rng('default');
        P = 16;
        N = 1000;
        
        o.Z = BeBP(N,P/2,1)';
        o.K = size(o.Z,1);
        
        o.D = 0.5 + randn(P,o.K)/2 ;
        
        o.gs = 2;
        o.S = o.gs.^-0.5*randn(size(o.Z));

        o.X = o.D*(o.S.*o.Z);
        
        eps = 0.05;
        o.ge = 1/eps;
        o.Y = o.X + eps^0.5*randn(P,N);
    end
end

methods(TestMethodTeardown)
    function teardown(o)
        warning('on', 'MATLAB:structOnObject');
    end
end

methods (Test)
    function infer_ge(o)
        args = struct(o);
        args.inferGe = true;
        bpfa = BPFA(args);
        bpfa.learn();

        for i = 1:100, ge(i) = bpfa.sample_ge(); end

        o.verifyEqual(mean(ge), o.ge, 'RelTol', 5e-2);
    end

    function infer_gs(o)
        args = struct(o);
        args.inferGs = true;
        args.gs = 5; %update depends on current, adjusted away from actual
        bpfa = BPFA(args);
        bpfa.learn();

        for i = 1:100, gs(i) = bpfa.sample_gs(); end

        o.verifyEqual(mean(gs), o.gs, 'RelTol', 1e-2);
    end

    function infer_pi(o)
        args = struct(o);
        bpfa = BPFA(args);
        bpfa.learn();

        for i = 1:100, pie(:,i) = bpfa.sample_pi(); end

        o.verifyEqual(mean(pie,2), mean(o.Z,2), 'AbsTol', 1e-2);
    end

    function infer_D_noiseless(o)
        args = struct(o);
        rmfield(args,'D');
        args.inferD = true;
        args.inferGe = true;
        args.Y = o.X;
        args.ge = 1;
        bpfa = BPFA(args);
        bpfa.learn();

        err = o.X - bpfa.X;
        o.verifyEqual(sqrt(mean(err(:).^2)), 0, 'AbsTol', 1e-5);
    end

    function infer_D(o)
        args = struct(o);
        rmfield(args,'D');
        args.inferD = true;
        args.inferGe = true;
        args.ge = 1;
        bpfa = BPFA(args);
        bpfa.learn();

        err = o.X - bpfa.X;
        o.verifyEqual(sqrt(mean(err(:).^2)), 0, 'AbsTol', 1e-1);
    end
end

end

