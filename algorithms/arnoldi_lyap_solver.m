% Input : A E where A^T P E + E^T P A < 0
%         R initial matrix
% Output: P
% Calculating P matrix by Arnoldi algorithm (AQ=QH)
function [P] = arnoldi_lyap_solver(A,E,R,q)
    A_bar=inv(A)*E;

    [Q,H]=restarted_arnoldi(A_bar',R,q);
    [V,D]=eig(full(H));
    D=diag(D);
    while(~all(real(D)<=0))
        disp('postive poles detected, restart arnoldi procedures')
        [Q,H]=restarted_arnoldi(A_bar',R,size(H,1));
        [V,D]=eig(full(H));
        D=diag(D);
    end

    % [Q,H]=block_arnoldi(A_bar',R,q);
    % [V,D]=eig(full(H));

    disp('The rank of P matrix is');
    disp(numel(D));
    keyboard;
    P = inv(A)'*Q*V*V'*Q'*inv(A);
    assert(norm(real(P))/norm(imag(P))>1e10);
    P = real(P);
%
