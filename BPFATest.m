classdef BPFATest < matlab.unittest.TestCase

properties
    Y,X
    D,S,Z
    ge,gs

    K

    inferA = false
    inferD = false, inferS = false, inferZ = false
    inferGe = false,inferGs = false
end

methods (TestMethodSetup)
    function setup(o)
        rng('default');
        P = 64;
        N = 10000;
        
        o.Z = BeBP(N,P/2,1)';
        o.K = size(o.Z,1);
        
        o.D = binornd(255,0.5,P,o.K)/255;
        
        o.gs = 1;
        o.S = o.gs*randn(size(o.Z));

        o.X = o.D*(o.S.*o.Z);
        
        eps = 0.05;
        o.ge = 1/eps;
        o.Y = o.X + eps^0.5*randn(P,N);
    end
end

%methods(TestMethodTeardown)
%    function teardown(o)
%        
%    end
%end

methods (Test)
    function infer_ge(o)
        args = struct(o);
        args.inferGe = true;
        bpfa = BPFA(args);
        bpfa.learn();

        for i = 1:100, ge(i) = bpfa.sample_ge(); end

        o.verifyEqual(mean(ge), o.ge, 'AbsTol', 1e-1);
    end

    function infer_gs(o)
        args = struct(o);
        args.inferGs = true;
        args.gs = 5;
        bpfa = BPFA(args);
        bpfa.learn();

        for i = 1:100, gs(i) = bpfa.sample_gs(); end

        o.verifyEqual(mean(gs), o.gs, 'AbsTol', 5e-2);
    end

    function infer_pi(o)
        args = struct(o);
        bpfa = BPFA(args);
        bpfa.learn();

        for i = 1:100, pie(:,i) = bpfa.sample_pi(); end

        o.verifyEqual(mean(pie,2), mean(o.Z,2), 'AbsTol', 1e-2);
    end
end

end

