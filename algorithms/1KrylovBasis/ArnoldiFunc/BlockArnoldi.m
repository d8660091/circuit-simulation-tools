function [ Q, H] = BlockArnoldi( G, C, R, q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function:  Block Arnoldi Process with Re-Orthogonalization
%
% Block Arnoldi Algorithm for (partial) reduction to Hessenberg form
%       Determining the orthogonal basis Q of the Krylov subspace
%       Kr(A,R,n) where n < m (# of rows) .
%
% Inputs:
%     q: The number of the column vectords in Q matrix (or the Order of the
%        reduced system)
%      Q_{n*q} = kr(A, R, K) where A_{n*n} and B_{n*N}
%      n: size of the MNA
%      N: Number of the ports
%      K=ceil(q/N) or
%      accurately: q=k*N+l, here k=floor(q/N) and l=q-kN   see proceeding
%      paper page 722
%
% N.B.:  The reduced system of order q preserves the first K=floor(q/N)
% number of the block moments of the original net; where N is the number of
% the ports.
% This implies that for a desired predefined accuracy, the order of the
% reduced system should be increased with the increase in the number of
% ports.
%
% {NOTE1:} The major limitation of the block Arnoldi is the
% orthogonality loss that occurs between the Krylov vectors as J (index of
% the order) increases [1, PP 219]
% For the record I have practically faced this many times even in a n-piseg
% of RC netwrok when higher order is required.
%
% Reference:
% This is an implementation of the Pseudocode of the BLOCK Arnoldi Based
%   PRIMA Algorithm outlined in:
% [1] M. Celick, et. al.,"IC Interconnect Analysis", 2002, KluwerAcademeic
% Publishers, PP 219
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE that Q is expected to be a dense matrix, then using sparse Arithmetic
%(by sparse R and A) to work out Q will make this Algorithm 200 times slower;
% however the the result would be more accurate rather than the one from
% defining matrices R and A as shown below:

% R=full(R); %to speed up
% A=full(A); %to speed up

[n, N] = size(R);  %n: MNAsize & N: PortsNr
NrBlkM = ceil(q/N);  %About to form kr(A, R, NrBlkM)
q = NrBlkM*N;

H = zeros(q,q);
Q = zeros(n,q);

[Q(:,1:N),~] = qr(full(R) ,0);

for j = 1: NrBlkM
%     fprintf('\n Block-Moment # %g',j);
    ci = [(j-1)*N+1 : j*N];  %column of the location the Block in H
    Z = -G\(C*Q(:,ci));
    for eps = 1:2     %Double Orthogonalization {see the "NOTE1" above}

        %%% Modified Gram-Schmidt orthogonalization:
        for i = 1:j   %Moves along the row in H
          ri = [(j-i)*N+1 : (j-i)*N+N ]; %row of the location the Block in H
          h = Q(:,ri).' * Z;
          Z = Z - Q(:,ri)* h;
          H(ri, ci) = H(ri, ci)+h;
        end

    end
    if i<NrBlkM
      [V,~] = qr(Z,0);
      rii = [1+j*N: (j+1)*N];
      Q(:, rii) = V;
      H(rii, ci)= h;
    end

end

end %EOF
