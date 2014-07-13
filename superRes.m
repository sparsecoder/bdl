function sr = superRes(img, scale)
[P,N] = size(img.patches);
srPatches = zeros(scale^2*P,N);
patches = img.patches;
for i = 1:N
    fprintf('%d/%d\n',i,N)
    [srPatch, srPatch0, lambda] = LPL(reshape(patches(:,i), img.n1, img.n2), scale);
    srPatches(:,i) = srPatch(:);
end

srImg = zeros(scale*size(img.img));
srImg = Image(srImg, struct('n1',scale*img.n1,'n2',scale*img.n2,'d1',scale*img.d1,'d2',scale*img.d2));

sr = srImg.recon(srPatches);

