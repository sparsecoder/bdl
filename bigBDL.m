function x = bigBDL(D,y, l1,l2)
x0 = D\y;
res = 1e8;
dff = 1e8;
it = 1;
while dff > 1e-4 %|| it < 20; 
    xOld = x;

    x = x0./(l2 + l1./x);
    %xH = xH0.*xH./(lambda + xH);
    %xH = (scale/2)*xH0 - xH0./(1 + xH/lambda)

    resOld = res;
    %if(mod(it,1) == 0)
        %FxH = fft2(xH);
        res = norm(y - D*x)
        res = sqrt(mean((xOld(:) - x(:)).^2))
        dff = resOld - res;
    %end

    it = it + 1
end