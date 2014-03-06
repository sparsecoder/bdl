P = 32;

D = normalize([dct([eye(P) 0.05*randn(P)]) eye(P)]);
K = size(D,2);
Z = BeBP(K,100,1);
N = size(Z,2);
%Z = ones(K,N);
S = randn(size(Z));
X = D*(S.*Z) + 0*randn(P,N);

img = Image('texture.jpg');

bpfa = BPFA(img.patches, 256);

% D works
%bpfa.D = D;
%bpfa.sampleD = false;

% S kind of works
%bpfa.S = S;
%bpfa.sampleS = false;

% Z doesnt really work
%bpfa.Z = Z;
%bpfa.sampleZ = false;
%bpfa.sampleA = false;
    
bpfa.learn(2500);
