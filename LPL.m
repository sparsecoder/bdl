function xH = LPL(xL, scale, lambda)

[n,~] = size(xL); %square...
N = round(n*scale);

Z = eye(N);
Z(n/2+1 : N-n/2, :) = 0;
W = dftmtx(N);
DTD =  n/N^2*(W'*Z*W);

FxL = fft2(xL);
FxH0 = zeros(N,N);
IJ = [1:n/2, N-n/2+1:N];
FxH0(IJ,IJ) = FxL;
xH0 = abs(ifft2(FxH0));
xH0 = xH0 - min(xH0(:));
xH0 = xH0/max(xH0(:));

xH = xH0;
lambda = .25*ones(N,N);
imagesc(abs(reshape(xH,N,N))); title(0); colormap('gray'); drawnow;
for it = 1:4
    xH = xH - min(xH(:));
    xH = xH/max(xH(:));
    %     xH = 0.25*xH0 + xH;
    %     imadjust(xH);
    
    imagesc(xH); title(it); drawnow;
    
    xH = blockInv(xH,DTD,xH0,lambda,N);
    xH = blockInv(xH',DTD',xH0',lambda',N)';
    
end
end

function xH = blockInv(xH,DTD,xH0,lambda,N)
parfor i=1:N
    tic
    %         x = diag(lambda(:,i).*xH(:,i).^-1);
    %         xH(:,i) = abs(pcg(DTD + x, xH0(:,i), [],[], x));
    
    
    x = xH(:,i);
    x2 = sqrt(x);
    A = diag(lambda(:,i)) + bsxfun(@times,bsxfun(@times,DTD,x2),x2');
    xH(:,i) = abs(x2.*pcg(A, x2.*xH0(:,i), [],[], diag(1 + x)));
    
    fprintf('%d\t%f\n',i, toc)
end
end
