P = 64;
N = 10000;

Z = BeBP(N,P/2,1)';
K = size(Z,1);
D = normalize(binornd(255,0.5,P,K)/255);
K = size(D,2);
%Z = ones(K,N);
S = randn(size(Z));
X = D*(S.*Z);
eps = 0.05;
X = normalize(X);
Y = X + eps^0.5*randn(P,N);

img = Image('barbara256.png');

%bpfa = BPFA(Y, X, K);
bpfa = BPFA(img.patches, img.patches0, 128)

%bpfa.D = D; bpfa.sampleD = false;
%bpfa.ge = 1/eps;

%bpfa.S = S; bpfa.sampleS = false;
%bpfa.gs = 1;
%bpfa.Z = Z; bpfa.sampleZ = false;
%bpfa.pie = sum(Z,2)/N;

%bpfa.init(bpfa)   
bpfa.learn(100);

