path = 'stem/5nm/';
files = dir([path '*.tif']);
nPatch = 1000;
psize = 8^2;
p.n1 = psize^0.5;
p.n2 = psize^0.5;
p.random = nPatch;
patches = zeros(psize, nPatch*length(files));
for i=1:length(files)
    fprintf('loading %s\n', files(i).name);
    img = Image([path files(i).name], p);
    patches(:,(i-1)*nPatch+1:i*nPatch) = img.patches;
end

b = BPFA(patches, patches, 324)
b.learn(100);

