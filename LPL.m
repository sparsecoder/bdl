function [xH, xH0, lambda] = LPL(xL, scale)

    [n,~] = size(xL); %square...
    N = round(n*scale);
    
    FxL = fft2(xL);
    FxLnorm = norm(FxL(:));
    FxH0 = zeros(N,N);
    IJ = [1:n/2, N-n/2+1:N];
    FxH0(IJ,IJ) = FxL;
    xH0 = abs(ifft2(FxH0));
    xH0 = xH0 - min(xH0(:));
    xH0 = xH0/max(xH0(:));
    xH0 = max(1e-6, min(1-1e-6, xH0));

    xH = xH0;
    lambda = 0.2;

    imagesc(abs(reshape(xH,N,N))); title(0); colormap('gray'); drawnow;
    res = 1e8;
    dff = 1e8;
    it = 1;
    while dff > 1e-4 %|| it < 20; 
        xHold = xH;

        xH = xH0./(scale + lambda./xH);
        %xH = xH0.*xH./(lambda + xH);
        %xH = (scale/2)*xH0 - xH0./(1 + xH/lambda)
        
        resOld = res;
        %if(mod(it,1) == 0)
            %FxH = fft2(xH);
            %res = norm(FxL - FxH(IJ,IJ), 'fro')/FxLnorm %+ sum(abs(xH(:)))
            res = sqrt(mean((xHold(:) - xH(:)).^2))
            dff = resOld - res
        %end

        imagesc(abs(xH)); title(sprintf('%d: %f',it, dff)); drawnow;pause(0.5)
        it = it + 1
    end
    xH = xHold;
    %imagesc(abs(xH)); drawnow;
end
