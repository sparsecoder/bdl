P = 7;
K = 4;
N = 50;

D = normalize(randn(P,K));
Z = BeBP(K,1500,1);
Z = Z(:,1:N);
S = normalize(randn(size(Z)));

X = D*(S.*Z) + 10e-4*randn(P,N);

bpfa = BPFA(X,K);
% bpfa.D = D;
% bpfa.sampleD = false;
bpfa.S = S;
bpfa.Z = Z;
bpfa.sampleA = false;

bpfa.learn(300);
