% Input : os = original system
%         U  = left projection matrix
%         V  = right projection matrix
%         q  = order of reduced system
% Output:
function [rs] = reduce_arnoldi(os,q);
    AE=inv(os.A)*os.E;
    V=block_arnoldi(AE,-inv(os.A)*os.B,q);
    rs.E=V'*os.E*V;
    rs.A=V'*os.A*V;
    rs.B=V'*os.B;
    rs.C=os.C*V;
%
