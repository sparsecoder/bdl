%lr = im2double(imread('stem/BSTO-Al2O3-62-PV_5_20.tif'));
lr = im2double(imread('stem/AlN-GaN-NH3_1_20.tif'));
lr = lr(1:1024,:);

lrRegion = im2double(imread('stem/AlN-GaN-NH3_1_20ROI.tif'));
lrRegion = lrRegion(15:end-15, 15:end-15);

%hr = im2double(imread('stem/BSTO-Al2O3-62-PV_7_5.tif'));
hr = im2double(imread('stem/AlN-GaN-NH3_6_5.tif'));
hr = hr(1:1024,:);

scale = 1024/size(lrRegion,1);

% [lrROI, hrAlign, tform] = alignImages(hr.^2,lr.^1.5, true);
% 
% hrAlign = sqrt(hrAlign);
% lrRegion = lr(lrROI(1,2):lrROI(2,2),lrROI(1,1):lrROI(2,1));


%bicubic
lrNN = imresize(lrRegion, scale, 'nearest');
lrBC = imresize(lrRegion, scale, 'bicubic');


patchSize = 4;
pL.n1 = patchSize; pL.n2 = pL.n1;
pH.n1 = round(pL.n1*scale); pH.n2 = round(pL.n2*scale);
pL.d1 = 1; pL.d2 = pL.d1;
pH.d1 = 4; pH.d2 = pH.d1;
hrTrain = Image(hr,pH);
lrTrain = Image(lrRegion,pL);

rL = pL;
rH = pH;
rL.d1 = 2; rL.d2 = rL.d1;
rH.d1 = 8; rH.d2 = rH.d1;
hrTrain2 = Image(hr,rH);
lrTrain2 = Image(lrRegion,rL);

sL = pL;
sH = pH;
sL.d1 = 4; sL.d2 = sL.d1;
sH.d1 = 16; sH.d2 = sH.d1;
hrTrain3 = Image(hr,sH);
lrTrain3 = Image(lrRegion,sL);

hrTest = Image(zeros(round(scale*size(lr))), rH);
lrTest = Image(lr, rL);

K = 512;
iter = 20;
[lP,~] = size(lrTrain.patches);
[hP,N] = size(hrTrain.patches);

%stacked
patches = [lrTrain.patches; zeros(hP,N)] + [zeros(lP,N); hrTrain.patches];
[D, S, Z, ge, gs, Pi] = foo(patches,K, iter);
DL = D(1:lP,:);
DH = D(lP+1:end,:);

[~, Strain, Ztrain, ~, ~, ~] = foo(lrTrain2.patches,K, iter, DL, [], gs, Pi);
lrStkTrain = hrTrain2.recon(DH*(Strain.*Ztrain));
imwrite(lrStkTrain, 'lrStkTrain.pdf', 'pdf')

[~, Stest, Ztest, ~, ~, ~] = foo(lrTest.patches,K, iter, DL, [], gs, Pi);
lrStkTest = hrTest.recon(DH*(Stest.*Ztest));
imwrite(lrStkTest, 'lrStkTest.pdf', 'pdf')

%regressed
[DHR, SH, ZH, geRH, gsRH, PiRH] = foo(hrTrain.patches,K, iter);
[DLR, SL, ZL, geRL, gsRL, PiRL] = foo(lrTrain.patches,K, iter);
W = (SH.*ZH)/(SL.*ZL);
W2 = (SL.*ZL)/(SH.*ZH);

[~, StrainR, ZtrainR, ~, ~, ~] = foo(lrTrain2.patches,K, iter, DLR*W2, [], gsRH, PiRH);
lrRegTrain = hrTrain2.recon( DHR*(StrainR.*ZtrainR) );
imshow(lrRegTrain, 'Border', 'tight');
saveas(gcf,'lrReg.pdf')

imwrite(lrRegTrain, 'lrRegTrain.pdf', 'pdf')

[~, StestR, ZtestR, ~, ~, ~] = foo(lrTest.patches,K, iter, false, DLR*W2, [], gsR, PiR);
lrRegTest = hrTest.recon( DHR*(StestR.*ZtestR) );
imwrite(lrRegTest, 'lrRegTest.pdf', 'pdf')


% [D0, S0, Z0] = BPFA_GP(lrI.patches,ones(size(lrI.patches)),K,30);
% imagesc(lrI.recon(D0*(S0.*Z0)'))

% [D, S, Z] = BPFA_GP(patches,ones(size(patches)),K,30);
% A = (S.*Z)';
% Dl = D(1:lP,:);
% Dh = D(lP+1:end,:);
% 
% sqrt(sum(sum((lrI.patches - Dl*A).^2))/numel(A))
% sqrt(sum(sum((hrI.patches - Dh*A).^2))/numel(A))
% sqrt(sum(sum((lrI.img - lrI.recon(Dl*A)).^2))/numel(A))
% sqrt(sum(sum((hrI.img - hrI.recon(Dh*A)).^2))/numel(A))
% 
% [~, Sl, Zl] = BPFA_GP(lrI.patches,ones(size(lrI.patches)),K,30, false, Dl);
% Al = (Sl.*Zl)';
% sqrt(sum(sum((lrI.patches - Dl*Al).^2))/numel(A))
% sqrt(sum(sum((hrI.patches - Dh*Al).^2))/numel(A))
% sqrt(sum(sum((lrI.img - lrI.recon(Dl*Al)).^2))/numel(A))
% sqrt(sum(sum((hrI.img - hrI.recon(Dh*Al)).^2))/numel(A))
% 
% 
% [Dl, Sl, Zl] = BPFA_GP(lrI.patches,ones(size(lrI.patches)),K,30);
% [Dh, Sh, Zh] = BPFA_GP(hrI.patches,ones(size(hrI.patches)),K,30);
% % [~, WS, WZ] = BPFA_GP(Sl.*Zl,ones(size(Sl.*Zl)),K,300, false, Sh.*Zh);
% W = (Sh.*Zh)'/(Sl.*Zl)';
% imagesc(hrI.recon(Dh*W*(Sl.*Zl)'))
% 
% lrTest = Image(lr(300:700,624:1024), pL);
% hrTest = Image(zeros(size(lrTest.img) + pH.n1-pL.n1), pH);
% [~, Stest, Ztest] = BPFA_GP(lrI.patches,ones(size(lrI.patches)),K,30, false, Dl);
% Atest = (Stest.*Ztest)';
% imagesc(hrTest.recon(Dh*W*Atest))
% 
% 
% 
% pL.n1 = 4; pL.n2 = pL.n1;
% pH.n1 = round(pL.n1/tform.scale); pH.n2 = round(pL.n2/tform.scale);
% pL.d1 = 1; pL.d2 = pL.d1;
% pH.d1 = pL.d1; pH.d2 = pL.d2;
% hrTrain = Image(hrAlign,pH);
% lrTrain = Image(lrRegion(1:end-(pH.n1-pL.n1), 1:end-(pH.n2-pL.n2)),pL);
% 
% pL.n1 = 4; pL.n2 = pL.n1;
% pH.n1 = round(pL.n1/tform.scale); pH.n2 = round(pL.n2/tform.scale);
% pL.d1 = 4; pL.d2 = pL.d1;
% pH.d1 = pL.d1; pH.d2 = pL.d2;
% hrI = Image(hrAlign,pH);
% lrI = Image(lrRegion(1:end-(pH.n1-pL.n1), 1:end-(pH.n2-pL.n2)),pL);
% 
% K = 512;
% [D,DL, ~,~] = BPFA_SR(hrTrain.patches, lrTrain.patches,K,30);
% W = lrTrain.patches/hrTrain.patches;
% [~,~, S, Z] = BPFA_SR(lrI.patches, [] ,K,30, false, W*D);
% A = S'.*Z';
% figure(1); imagesc(hrI.recon(D*A))
% 
% [~,~, S, Z] = BPFA_SR(lrTrain.patches, [] ,K,30, false, DL);
% A = S'.*Z';
% figure(2); imagesc(hrTrain.recon(D*A))
% 
% 
% 
% %{
% Ah = BPFA(hrPatches, hrPatches, K);
% Ah.D = D.D;
% Ah.init(Ah);
% Ah.sampleD = false;
% Ah.learn(20);
% 
% alphaMap = BPFA(Ah.S.*Ah.Z, Ah.S.*Ah.Z, K);
% alphaMap.S = Al.S.*Al.Z;
% alphaMap.Z = alphaMap.S;
% alphaMap.Z(abs(alphaMap.Z)>0) = 1;
% alphaMap.sampleA = false;
% alphaMap.init(alphaMap);
% 
% alphaMap.learn(50);
% 
% inferHR = D.D*((Al.S.*Al.Z));
% 
% hrD.learn(20);
% lrD.learn(20);
% 
% alphaMap = BPFA(hrD.S.*hrD.Z, hrD.S.*hrD.Z, K);
% alphaMap.S = lrD.S.*lrD.Z;
% alphaMap.Z = alphaMap.S;
% alphaMap.Z(abs(alphaMap.Z)>0) = 1;
% alphaMap.sampleA = false;
% alphaMap.init(alphaMap);
% 
% alphaMap.learn(50);
% 
% inferHR = hrD.D*(alphaMap.D*lrD.A);
% imagesc(hrI.recon(inferHR))
%}