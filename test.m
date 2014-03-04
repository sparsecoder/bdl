P = 16;
K = 128;
N = 1000;

D = normalize(randn(P,K));
Z = BeBP(K,N,1);
Z = Z(:,1:N);
%Z = ones(K,N);
S = normalize(randn(size(Z)));
X = D*(S.*Z) + 1e-4*randn(P,N);

img = Image(dct2(eye(200)) + 1e-4*randn(200,200));


bpfa = BPFA(X,K);

% D works
%bpfa.D = D;
%bpfa.sampleD = false;

% S kind of works
bpfa.S = S;
bpfa.sampleS = false;

% Z doesnt really work
bpfa.Z = Z;
%bpfa.sampleZ = false;
bpfa.sampleA = false;
    
bpfa.learn(20);
