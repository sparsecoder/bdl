classdef Image<handle
properties
    img
    N
    n1,n2,d1,d2
end

methods
    function o = Image(img, param)
        if isnumeric(img), o.img = img; else o.img = imread(img); end
        o.img = im2double(o.img);
        o.N = size(o.img);
        if nargin>1
            if all(isfield(param, {'n1','n2'})), o.n1 = param.n1; o.n2 = param.n2;
            else o.n1 = 8; o.n2 = 8; end
            if all(isfield(param, {'d1','d2'})), o.d1 = param.d1; o.d2 = param.d2;
            else o.d1 = 1; o.d2 = 1; end
        else
            o.n1 = 8; o.n2 = 8; o.d1 = 1; o.d2 = 1;
        end
    end
    
    
    
    function patches = patch(o)
        [I,J] = o.getCoords();
        K = length(I);
        
        patches = zeros(o.n1*o.n2*o.N(3), K);
        for k=1:K
            p = o.img(I(k):I(k)+o.n1-1, J(k):J(k)+o.n2-1, :);
            patches(:,k) = p(:);
        end
    end
    
    function img = recon(o, patches)
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