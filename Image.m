classdef Image<handle
properties
    img
    img0
    N
    n1,n2,d1,d2
    Mu
    Sigma
    patches
    patches0
end

methods
    function o = Image(img, param)
        if isnumeric(img), o.img = img; else o.img = imread(img); end

        o.img = rgb2gray(o.img);
        o.img = im2double(o.img);
        o.Mu = mean(o.img(:));
        o.img = o.img - o.Mu;
        o.Sigma = std(o.img(:));
        o.img = o.img/o.Sigma;
        %o.img = imresize(o.img, [nan, 300]);
        o.img0 = o.img;

        o.N = size(o.img);
        o.img = max(0,min(1, o.img0 + 0.2*randn(o.N) ));

        if length(o.N)==2, o.N(3) = 1; end
        if nargin>1
            if all(isfield(param, {'n1','n2'})), o.n1 = param.n1; o.n2 = param.n2;
            else o.n1 = 8; o.n2 = 8; end
            if all(isfield(param, {'d1','d2'})), o.d1 = param.d1; o.d2 = param.d2;
            else o.d1 = 1; o.d2 = 1; end
        else
            o.n1 = 8; o.n2 = 8; o.d1 = 2; o.d2 = 2;
        end
        
        o.patches = o.patch(o.img);
        o.patches0 = o.patch(o.img0);
        % o.normalize();
    end
    
    function patches = patch(o, x)
        [I,J] = o.getCoords();
        K = length(I);
        
        patches = zeros(o.n1*o.n2*o.N(3), K);
        for k=1:K
            p = x(I(k):I(k)+o.n1-1, J(k):J(k)+o.n2-1, :);
            patches(:,k) = p(:);
        end
    end
    
    function img = recon(o, patches, lambda)
        [I,J] = o.getCoords();
        img = zeros(o.N);
        w = zeros(o.N);
        K = length(I);
        
        for k = 1 : K
            i = I(k):I(k)+o.n1-1;
            j = J(k):J(k)+o.n2-1;
            img(i, j,:) = img(i, j,:) + reshape(patches(:, k), [o.n1, o.n2, o.N(3)]);
            w(i, j,:) = w(i, j,:) + 1;
        end

        img = img ./ w;
        if nargin>2
            img = img + lambda*o.img;
        end
        img = img*o.Sigma + o.Mu;
    end
    
    function [I,J] = getCoords(o,I,J)
        if nargin < 3
            I = uint16(1 : o.d1 : o.N(1) - o.n1 + 1);
            J = uint16(1 : o.d2 : o.N(2) - o.n2 + 1);
        end
        [J, I] = meshgrid(J, I);
        J = J(:);
        I = I(:);
    end
    
end
end
