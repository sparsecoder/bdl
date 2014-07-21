function [xH, xH0, lambda] = LPL(xL, scale)

    [n,~] = size(xL); %square...
    N = round(n*scale);

Zdown = @(x) x([1:n/2, N-n/2+1:N], [1:n/2, N-n/2+1:N]);
Fn = @(x)reshape(Zdown(fft2(reshape(x,N,N))),n*n,1);
F0 = sparse(N,N);
Zup = @(x) [x(1:n/2,1:n/2),     zeros(n/2, N-n),  x(1:n/2,n/2+1:end);
                                zeros(N-n, N);
            x(n/2+1:end,1:n/2)  zeros(n/2, N-n)   x(1:n/2,n/2+1:end) ];
Fnt = @(xhat)reshape(ifft2(Zup(reshape(xhat,n,n))),N*N,1);

    Z = speye(N);
    Z(n/2+1 : N-n/2, :) = 0;
    Z = kron(Z,Z);
    Z = reshape(diag(Z),N,N);

    FxL = fft2(xL);
    FxH0 = zeros(N,N);
    IJ = [1:n/2, N-n/2+1:N];
    FxH0(IJ,IJ) = FxL;
    xH0 = abs(ifft2(FxH0));
    xH0 = xH0 - min(xH0(:));
    xH0 = xH0/max(xH0(:));
    xH0 = max(1e-6, min(1-1e-6, xH0));

    xH = xH0;%randn(size(xH0));
    lambda = 0.2;%*betapdf(xH0,0.1,0.1);%1/lambda0,lambda0);

    %imagesc(abs(reshape(xH,N,N))); title(0); colormap('gray'); drawnow;
    res = 1e8;
    dff = 1e8;
    it = 1;
    while dff > 1e-1% || it < 50; 
        xHold = xH;

        %x = xH(:);
        %x2 = sqrt(x);
        %DtDX2 = Fnt(Fn(x2));
        %A = lambda(:) + x2.*DtDX2; %X^.5*DTD*X^.5
       % [xH(:,i),~] = pcg(A, x2.*xH0(:), [],[], diag(1+x));
%        xH = reshape(x.*xH0(:)./A, N,N);
        %xH = xH0./(lambda + lambda./xH);
       xH = xH0.*xH./(lambda + xH);
        
        %if mod(it,2)
        %    xH = blockInv(xH,DTD,xH0,lambda,N);
        %else
        %    xH = blockInv(xH',DTD,xH0',lambda',N)';
        %end

        resOld = res;
        if(mod(it,100) == 0)
            FxH = fft2(xH);
            res = norm(FxL - FxH(IJ,IJ), 'fro')
            dff = resOld - res
        end

        %imagesc(abs(xH)); title(sprintf('%d: %f',it, diff)); drawnow;
        it = it + 1
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
