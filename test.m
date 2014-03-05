P = 16;

D = normalize([dct([eye(P) 0.05*randn(P)]) eye(P)]);
K = size(D,2);
Z = BeBP(K,100,1);
N = size(Z,2);
%Z = ones(K,N);
S = randn(size(Z));
X = D*(S.*Z) + 1e-2*randn(P,N);


%img = imread('skyline.jpg');
%img = 0.2989*img(:,:,1) + 0.5870*img(:,:,2) + 0.1140*img(:,:,3);
%img = Image(img);

bpfa = BPFA(X, K);

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
    
bpfa.learn(100);
