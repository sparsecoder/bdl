lr = im2double(imread('stem/BSTO-Al2O3-62-PV_5_20.tif'));
lr = lr(1:1024,:);

hr = im2double(imread('stem/BSTO-Al2O3-62-PV_7_5.tif'));
hr = hr(1:1024,:);

[lrROI, hrAlign, tform] = alignImages(hr.^2,lr.^1.5, true);

hrAlign = sqrt(hrAlign);
lrRegion = lr(lrROI(1,1):lrROI(2,1), lrROI(1,2):lrROI(2,2));

pL.n1 = 1; pL.n2 = 1; pL.d1 = 1; pL.d2 = 1;
pH.n1 = round(pL.n1/tform.scale); pH.n2 = round(pL.n2/tform.scale);
pH.d1 = 1; pH.d2 = 1;
hrI = Image(hrAlign,pH);
lrI = Image(lrRegion(1:end-(pH.n1-pL.n1), 1:end-(pH.n2-pL.n2)),pL);

K = 512;
[lP,N] = size(lrI.patches);
[hP,N] = size(hrI.patches);
lrPatches = [lrI.patches; zeros(hP,N)];
hrPatches = [zeros(lP,N); hrI.patches];
patches = [lrPatches + hrPatches];
D = BPFA(patches, patches, K);
D.learn(20);
Dl = D.D(1:lP,:);
Dh = D.D(lP+1:end,:);

norm(lrI.patches - Dl*D.A)
norm(hrI.patches - Dh*D.A)

Al = BPFA(lrI.patches, lrI.patches, K);
Al.D = Dl;
Al.init(Al);
Al.sampleD = false;
Al.learn(20);

norm(lrI.patches - Dl*Al.A)
norm(hrI.patches - Dh*Al.A)


%{
Ah = BPFA(hrPatches, hrPatches, K);
Ah.D = D.D;
Ah.init(Ah);
Ah.sampleD = false;
Ah.learn(20);

alphaMap = BPFA(Ah.S.*Ah.Z, Ah.S.*Ah.Z, K);
alphaMap.S = Al.S.*Al.Z;
alphaMap.Z = alphaMap.S;
alphaMap.Z(abs(alphaMap.Z)>0) = 1;
alphaMap.sampleA = false;
alphaMap.init(alphaMap);

alphaMap.learn(50);

inferHR = D.D*((Al.S.*Al.Z));

hrD.learn(20);
lrD.learn(20);

alphaMap = BPFA(hrD.S.*hrD.Z, hrD.S.*hrD.Z, K);
alphaMap.S = lrD.S.*lrD.Z;
alphaMap.Z = alphaMap.S;
alphaMap.Z(abs(alphaMap.Z)>0) = 1;
alphaMap.sampleA = false;
alphaMap.init(alphaMap);

alphaMap.learn(50);

inferHR = hrD.D*(alphaMap.D*lrD.A);
imagesc(hrI.recon(inferHR))
%}