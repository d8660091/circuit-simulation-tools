% Input : os = original system
%         P  = Lyapunov solution
%         q  = order of reduced system
% Output:
function [rs] = reduce_stable(os,P,q);
    V=block_arnoldi(inv(os.A)*os.E,-inv(os.A)*os.B,q);
    U=(V'*os.E'*P)';

    U=orth(U);
    q=size(U,2);
    V=V(:,1:q);

    rs.E=U'*os.E*V;
    rs.A=U'*os.A*V;
    rs.B=U'*os.B;
    rs.C=os.C*V;

    % U=V'*os.E'*P;
    % rs.E=U*os.E*V;
    % rs.A=U*os.A*V;
    % rs.B=U*os.B;
    % rs.C=os.C*V;
