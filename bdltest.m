n = 4096;
x = binornd(1,0.1,n,1)*randn(n,1);
y = fft(x);
xhat = bigBDL(dftmtx(n), y, 0.2, 1);

norm(x - xhat);