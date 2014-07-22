addpath('NESTA_v1.1')
Setup_Nesta

fname = 'stem/GTO-LSAT-49_5_20.tif';
%fname = 'stem/GTO-LSAT-297_22_10.tif';
lr = im2double(imread(fname));
lr = lr(1:1024,:);

lr = lr((1:512)+256,(1:512));
%lr = lr(101:116+16, (101:116+16)-30);

% lrImg = Image(lr, struct('n1',16,'n2',16,'d1',8,'d2',8));

% sr = superRes(lrImg, 8);


scale = 10;
tic
[xH, xH0, lambda] = LPL(lr, scale);
toc
xH = abs(xH);
imagesc(xH)
%imwrite(xH/max(xH(:)), 'srLPL.png', 'png')


%[n,~] = size(lr); %square...
%N = round(n*scale);
%% Z = speye(N);imagesc
%% Z(n/2+1 : N-n/2,:) = [];
%% W = dftmtx(N)/n;
%% W = kron(Z*W,Z*W);
%
%Zdown = @(x) x([1:n/2, N-n/2+1:N], [1:n/2, N-n/2+1:N]);
%Fn = @(x)reshape(Zdown(fft2(reshape(x,N,N))),n*n,1);
%F0 = sparse(N,N);
%Zup = @(x) [x(1:n/2,1:n/2),     zeros(n/2, N-n),  x(1:n/2,n/2+1:end);
%                                zeros(N-n, N);
%            x(n/2+1:end,1:n/2)  zeros(n/2, N-n)   x(1:n/2,n/2+1:end) ];
%Fnt = @(xhat)reshape(ifft2(Zup(reshape(xhat,n,n))),N*N,1);
%
%opts = struct();
%opts.TOlVar = 1e-2;
%opts.verbose = 10;
%opts.maxintiter = 5;
%opts.maxiter = 400;
%opts.stoptest = 1;
%opts.typemin = 'tv';
%
%mu = 1e-6;
%tic
%yhat = fft2(lr);
%% x = NESTA(W,[],yhat(:),mu,eps,opts);
%[x,niter,residuals] = NESTA_UP(Fn,Fnt,yhat(:),0.2,n/2,mu,opts);
%toc
%X = abs(reshape(x,N,N));
%
%%figure(1)
%%imagesc(X)
%%figure(2)
%%imagesc(lr)
%
%imwrite(X/max(X(:)), 'sr.png', 'png')
