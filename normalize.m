function X = normalize(X)
    X = bsxfun(@minus, X, mean(X));
    X = bsxfun(@rdivide, X, std(X));
end