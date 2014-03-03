P = 120;
K = 20;
N = 50000;

%D = normalize(randn(P,K));
%Z = BeBP(K,N,1);
%Z = Z(:,1:N);
%S = normalize(randn(size(Z)));
%
%X = D*(S.*Z) + 1e-4*randn(P,N);

bpfa = BPFA(X,K);

% D works
%bpfa.D = D;
%bpfa.sampleD = false;

% S kind of works
%bpfa.S = S;
%bpfa.sampleS = false;

% Z doesnt really work
%bpfa.Z = Z;
%bpfa.sampleZ = false;


bpfa.learn(20);
