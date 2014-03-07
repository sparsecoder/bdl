P = 16;

D = normalize(binornd(255,0.5,P,64)/255);
K = size(D,2);
Z = BeBP(K,100,1);
N = size(Z,2);
%Z = ones(K,N);
S = randn(size(Z));
X = D*(S.*Z);
Y = X + 1e-1*randn(P,N);

img = Image('skyline.jpg');

bpfa = BPFA(Y, X, K);
%bpfa = BPFA(img.patches, img.patches0, 96);

%bpfa.D = D; bpfa.sampleD = false;

%bpfa.S = S; bpfa.sampleS = false;
bpfa.Z = Z; bpfa.sampleZ = false;

bpfa.init(bpfa)   
bpfa.learn(100);
