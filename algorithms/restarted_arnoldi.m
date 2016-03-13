% A*V = V*R + x,  all eigenvalues of H is less than 0
function [Q,H,ritz] = restarted_arnoldi(A,m,k);
    v=rand(size(A,1),1);
    n = length(v);
    H = zeros(k);
    V = zeros(n,k);

    v = v/norm(v);
    w = A*v;
    alpha = v'*w;

    V(:,1) = v; H(1,1) = alpha;
    f = w - v*alpha;

    for j = 2:k,
        beta = norm(f);
        v = f/beta;
        H(j,j-1) = beta;
        V(:,j)   = v;
        w = A*v;
        h = V(:,1:j)'*w;
        f = w - V(:,1:j)*h;
        H(1:j,j) = h;
    end
