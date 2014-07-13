function [xH, xH0, lambda] = LPL(xL, scale)

    [n,~] = size(xL); %square...
    N = round(n*scale);

    Z = eye(N);
    Z(n/2+1 : N-n/2, :) = 0;
    W = dftmtx(N);
    DTD =  1/N*(W'*Z*W);

    FxL = fft2(xL);
    FxH0 = zeros(N,N);
    IJ = [1:n/2, N-n/2+1:N];
    FxH0(IJ,IJ) = FxL;
    xH0 = abs(ifft2(FxH0));
    xH0 = xH0 - min(xH0(:));
    xH0 = xH0/max(xH0(:));
    xH0 = max(1e-6, min(1-1e-6, xH0));

    xH = xH0;
    lambda0 = 2;
    lambda = betapdf(xH0,1/lambda0,lambda0);

    %imagesc(abs(reshape(xH,N,N))); title(0); colormap('gray'); drawnow;
    del = 1;
    diff = 1;
    it=1;
    while diff > 0
        xHold = xH;

        if mod(it,2)
            xH = blockInv(xH,DTD,xH0,lambda,N);
        else
            xH = blockInv(xH',DTD,xH0',lambda',N)';
        end

        delOld = del;
        del = norm(xHold - xH)/norm(xH0);
        diff = delOld - del;

        %imagesc(abs(xH)); title(sprintf('%d: %f',it, diff)); drawnow;
        it = it + 1;
    end
    xH = xHold;
    %imagesc(abs(xH)); drawnow;
end

function xH = blockInv(xH,DTD,xH0,lambda,N)
    parfor i=1:N
        %         x = diag(lambda(:,i).*xH(:,i).^-1);
        %         xH(:,i) = abs(pcg(DTD + x, xH0(:,i), [],[], x));

        x = xH(:,i);
        x2 = sqrt(x);
        A = diag(lambda(:,i)) + bsxfun(@times,bsxfun(@times,  DTD,x2),x2'  ); %X^.5*DTD*X^.5
        [xH(:,i),~] = pcg(A, x2.*xH0(:,i), [],[], diag(1+x));
        xH(:,i) = x2.*xH(:,i);
    end
end
