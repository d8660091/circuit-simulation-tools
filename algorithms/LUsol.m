function X = LUsol(A,B)
%
% This is to solve linear equations in the form of: AX=BU using LU deomposition
%--------------------------------------------------------------------------
persistent sp;
if isempty(sp)
    sp=issparse(A);
end

if sp; 
   [L1, U1, P1, Q1 ]= lu( A, 0.1 );
   X = full( Q1*(U1 \( L1 \ (P1*B) )) );
else
   [L1, U1, P1, Q1 ]= lu( A );
   X = full(U1 \( L1 \ (P1*B) ) );
end
%EOF
