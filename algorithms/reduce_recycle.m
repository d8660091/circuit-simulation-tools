% Input : A where A^T P E + E^T P A < 0
%         q the size of reduced model
% Output: P
% Reduce the original stable system to a smaller stable system just use one Arnoldi process.
function [rs] = reduce_recycle(os,q)
    [Q,H]=restarted_arnoldi(os.A\os.E,-os.A\os.B,q);
    while(~all(real(eig(full(H)))<=0))
        disp('postive poles detected, restart arnoldi procedures')
        [Q,H]=restarted_arnoldi(os.A\os.E,-os.A\os.B,size(H,1));
    end
    [V,~]=eig(full(H'));
    P = inv(os.A)'*Q*V*V'*Q'*inv(os.A);
    V=Q; % original V will not be used.
    rs.E=V'*os.E'*P*os.E*V;
    rs.A=V'*os.E'*P*os.A*V;
    rs.B=V'*os.E'*P*os.B;
    rs.C=os.C*V;
%
