P = 7;
K = 4;
N = 500;

D = normalize(randn(P,K));
Z = BeBP(K,N,1);
Z = Z(:,1:N);
S = 255*normalize(randn(size(Z)));

X = D*(S.*Z) + 10e-6*randn(P,N);

bpfa = BPFA(X,K);

% D works
bpfa.D = D;
bpfa.sampleD = false;

% S kind of works
%bpfa.S = S;
%bpfa.sampleS = false;

% Z doesnt really work
bpfa.Z = Z;
bpfa.sampleZ = false;


bpfa.learn(100);
