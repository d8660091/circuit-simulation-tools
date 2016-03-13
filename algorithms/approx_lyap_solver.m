% Input : A E where A^T P E + E^T P A < 0
% Output: P
% Calculating P matrix by largest eigenvalues (AV=VD)
function [V,D_matrix,P] = approx_lyap_solver(os,q,p,B)
    A=os.A;
    E=os.E;
    B=os.B;

    % Another way
    % [Q,H]=block_arnoldi(A_bar',R,q);
    % [V,D]=eig(full(H));
    % P = inv(A)'*Q*V*V'*Q'*inv(A);

    % A_bar=A\E;
    % [V,D]=eigs(A_bar',q,'lm');

    [V,D]=eigs(E',A',q);
    %V=A'*V;

    % [V,D]=eigs(A',E',q);
    % V=A'*V;
    % D=1./D;
    D=diag(D);
    i=real(D)<0;
    D=D(i);
    V=V(:,i);

    % select largest real-part eigenvalues
    [~,i]=sort(real(D));
    i=i(1:p);
    D=D(i);
    V=V(:,i);
    % remove lonely complex eigenvectors again
    if(~isreal(D(end)) && D(end-1)~=conj(D(end)))
        V=V(:,1:p-1);
        D=D(1:p-1);
    end

    % random select eigenvectors
    % i=rand(numel(D),1);
    % i=i<0.5;
    % D=D(i);
    % V=V(:,i);
    %
    disp('The rank of P is');
    disp(numel(D));

    %% P=inv(A)'*V*V'*inv(A);
    P=0;
    if(nargout==3)
        tmp=A'\V;
        P=tmp*tmp';
        P=real(P);
    end
    D_matrix=diag(D);
