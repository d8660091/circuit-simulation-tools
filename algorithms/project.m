% Projecting original system to reduced system
% Input : os = original system
%         U  = left projection matrix
%         V  = right projection matrix
% Output:
function [rs] = reduce_stable(os,U,V);
    rs.E=U'*os.E*V;
    rs.A=U'*os.A*V;
    rs.B=U'*os.B;
    rs.C=os.C*V;
%
