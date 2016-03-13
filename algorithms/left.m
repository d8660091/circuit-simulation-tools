% Computes left projecting matrix V
% Input : os = original system
% Output: U = left projection matrix
%
function [U] = left(os,q);
    A_bar=inv(os.A)*os.E;
    % [Q,H]=block_arnoldi(A_bar',R,q);
    % [V,D]=eig(full(H));
    % P = inv(A)'*Q*V*V'*Q'*inv(A);

    [V,D]=eigs(A_bar',q);
    if(D(end-1)~=conj(D(end)))
        [V,D]=eigs(A_bar',q-1);
    end
    D=diag(D);
    i=real(D)<0;
    D=D(i);
    V=V(:,i);
    tmp=inv(os.A)'*V*V';
    U=orth(tmp);
