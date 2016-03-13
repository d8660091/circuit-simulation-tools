function [ Q, H] = block_arnoldi( A, R, q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function:  Block Arnoldi Process with Re-Orthogonalization
%
% Block Arnoldi Algorithm for (partial) reduction to Hessenberg form
%       Determining the orthogonal basis Q of the Krylov subspace
%       Kr(A,R,n) where n < m (# of rows) .
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, N] = size(R);  %n: MNAsize & N: PortsNr

NrKr = ceil(q/N);  %About to form kr(A, R, NrKr)

%H = zeros(q,q);
H = sparse(NrKr*N,NrKr*N);
[Q,X] = qr(R,0); clear X;


for j = 1: NrKr
    ci = [(j-1)*N+1 : j*N];  %column of the location the Block in H
    Z = A * Q(:,ci);         %<= using lu f/b is preferable here to improve cpu economy

    for eps = 1:2     %Double Orthogonalization

        %%% Modified Gram-Schmidt orthogonalization:
        for i = 1:j   %Moves along the row in H
          ri = [(j-i)*N+1 : (j-i)*N+N ]; %row of the location the Block in H
          h = Q(:,ri).' * Z;
          Z = Z - Q(:,ri)* h;
          H(ri, ci) = H(ri, ci)+h;
        end

    end
    if i<NrKr
      [V,h] = qr(Z,0);
      Q = [Q,V];
      ri = [1+j*N: (j+1)*N];
      H(ri, ci) = h;
    end
end
end %EOF
