function [V,T,f] = lanczos(A,v,k);
%
%   Input:  A -- an n by n matrix  (A = A' assumed)
%           v -- an n  vector (v .ne. 0 assumed)
%           k -- a positive integer (k << n assumed)
%
%   Output: V -- an n by k orthogonal matrix
%           T -- a k by k symmetric tridiagonal matrix
%           f -- an n vector
%
%
%           with   AV = VT + fe_k'
%
%   In real life you would not store V, and you would store T as two
%   vectors  a = diag(T) , b = diag(T,-1)
%
%   D.C. Sorensen
%   21 Feb 00
%
    n = length(v);
    T = zeros(k);
    V = zeros(n,k);

    v1 = v/norm(v);

    f = A*v1;
    alpha = v1'*f;
    f = f - v1*alpha;

    V(:,1) = v1; T(1,1) = alpha;

    for j = 2:k,

        beta = norm(f);
        v0 = v1; v1 = f/beta;

        f = A*v1 - v0*beta;
        alpha = v1'*f;
        f = f - v1*alpha;

        T(j,j-1) = beta; T(j-1,j) = beta; T(j,j) = alpha;
        V(:,j)   = v1;

    end



