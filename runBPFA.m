rng('default');

img = imread('barbara256.png');
img = double(img)/255;
img = img + 0.01*randn(size(img));
p = 6;
patches = im2col(img, [p,p], 'sliding');
K=92;
[D,S,Z,ge,gs] = foo(patches,K, 10);

%[D S Z] =BPFA_GP(patches, K, 4);