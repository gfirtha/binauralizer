function [K] = simplify_sphere(K,vertices)
eps = 1e-12;
for n = 1 :  1
    v0 = [vertices(1,K{n});vertices(2,K{n});vertices(3,K{n})];
    for m = size(v0,2) : -1 : 1
        R = sqrt(sum((v0(:,m) -v0(:,1:m)).^2,1));
        if any(R(1:m-1)<eps)
            K{n}(m) = [];
        end
    end
end
end

