%P = 64;
%K = 256;
%N = 100000;
%
%D = normalize(randn(P,K));
%Z = BeBP(K,N,1);
%Z = Z(:,1:N);
%%Z = ones(K,N);
%S = normalize(randn(size(Z)));
%X = D*(S.*Z) + 1e-4*randn(P,N);
%

img = Image('skyline.jpg');

bpfa = BPFA(img.patches, 512);

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
    
bpfa.learn(20);
