function [V,H] = qrstep(V,H,mu,k1,k2);
%
%   Input: V     -- a real square orthogonal matrix
%
%          H     -- a real square upper Hessenberg matrix
%
%          mu    -- a real  or complex shift
%
%          k1,k2 -- pointers to submatrix in k1:k2 block
%
%   Output: V    -- a real orthogonal matrix  
%
%           H    -- a real square upper Hessenberg matrix
%
%                   V <- VQ;   H <- Q'HQ;
%
%                   Q corresponds to a single real shift or 
%                     a double complex shift depending on mu
%
%   D.C. Sorensen
%   2 Mar 2000
%
    k = k2-k1+1;
    kr = k1:k2;
    n = length(V(:,1));
    
    eta = imag(mu);

    if (abs(eta) > 0), 
%                      mu is imaginary -- apply double shift
%
       xi = real(mu);
       [Q,R] = qr((H(kr,kr) - xi*eye(k))^2 + eta^2*eye(k));

    else  
%                      mu is real -- apply single shift
%
       [Q,R] = qr(H(kr,kr) - mu*eye(k));

    end

    H(kr,:) = Q'*H(kr,:);
    H(:,kr) = H(:,kr)*Q;
    V(:,kr) = V(:,kr)*Q;

%
%   clean up rounding error noise below first subdiagonal
%
    for j = k1:k2,
        H(j+2:n,j) = H(j+2:n,j)*0; 
    end
