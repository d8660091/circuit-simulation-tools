% Input : A
%         E
%         A^T P E + E^T P A < 0
% Output: P
function [P] = lyap_solver(A,E)
    A_bar=inv(A)*E;
    [V_bar,~]=eig(full(A_bar'));
    P = inv(A)'*V_bar*V_bar'*inv(A);
%
