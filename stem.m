files = dir('stem/*.tif');
nPatch = 10000;
patches = zeros(64, nPatch*length(files));
for i=1:length(files)
    img{i} = Image(['stem/' files(i).name]);
    inds{i} = randperm(size(img{i}.patches,2), 10000);
    patches(:,(i-1)*nPatch+1:i*nPatch) = img{i}.patches(:,inds{i});
end

b = BPFA(patches, patches, 512)
b.learn(100);

