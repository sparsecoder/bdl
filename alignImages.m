function [lrROI, hrAlign, tform] = alignImages(hr,lr, verbose)
if nargin < 3, verbose = false; end

[optimizer, metric]  = imregconfig('monomodal');
optimizer.MaximumIterations = 100;
p = floor(log2(min([size(hr) size(lr)]))) - 4;
[~,t] = imregister(hr,lr,'similarity', optimizer, metric, 'PyramidLevels', p);

theta = asin(t(1,2));
tform.theta = theta;
scale = 0.7*t(1,1)/cos(theta);
tform.scale = scale;
shift = t(3,1:2) + [75 45];

[M, N] = size(lr);
newt = [scale*cos(theta) sin(theta) 0;
        -sin(theta) scale*cos(theta) 0;
        shift 1];
tt = maketform('affine', newt);
reg = imtransform(hr, tt, 'XData', [1 N], 'YData', [1 M], 'Size', [M N]);
if verbose, imshowpair(lr,reg); end

crop = ceil(t(2,1)*1024) + 1;
tform.crop = crop;
hrRot = imrotate(hr, -theta*180/pi, 'bilinear');
hrRotCrop = hrRot(crop:end-crop, crop:end-crop);
hrAlign = imresize(hrRotCrop, scale);
cropShift = [shift(1)+crop/2 shift(2)];
tform.shift = cropShift;
newt = [scale 0 0;
        0 scale 0;
        cropShift 1];
tt = maketform('affine', newt);
reg = imtransform(hrRotCrop, tt, 'XData', [1 N], 'YData', [1 M], 'Size', [M N]);
if verbose, imshowpair(lr,reg); end

lrROI = [ceil(cropShift); ceil(cropShift) + size(hrAlign) - 1];
