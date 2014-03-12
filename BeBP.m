function [ Z ] = BeBP( N, a, b )
    c = poissrnd(a./(b+(0:N-1)));
    C = cumsum(c);
    ncol = C(end);
    Z = false(N,ncol);

    Z(1,1:c(1)) = 1;
    n = Z(1,:);
    for i=2:N
        Z(i, 1:C(i-1)) = rand(1,C(i-1)) < n(1:C(i-1))/(i + b-1);
        Z(i, C(i-1)+1:C(i)) = 1;
        n = n + Z(i,:);
    end
end

