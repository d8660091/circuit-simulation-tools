% Input: q the size of A
%        (q-r) the rank of E matrix
% Output: a random stable system
% Generating a random stable system
function os=random_sys(q,r)
    AE=sprand(q,q,0.03);
    [V,D]=eig(full(AE));
    % D=complex(-abs(real(D)),imag(D));

    tmp=diag(D);
    if(tmp(end-(r-1))==tmp(end-r))
        tmp(end-r:end)=0;
    else
        tmp(end-r+1:end)=0;
    end
    tmp=1e-4*tmp;
    % tmp(1:5)=1e-5*tmp(1:5);
    D=diag(tmp);

    AE=real(V*D*inv(V));
    A=sprand(q,q,0.03);
    E=A*AE;

    os.A=A;
    os.E=E;
    B=zeros(q,1);
    B(1)=1;
    os.B=B;
    os.C=B';
