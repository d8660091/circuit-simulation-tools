function R = FormR( Glu, B)
%
% Calculating A and R:
% R=(G^-1)B,  A=-(G^-1)C
%
% Recap:
% GX + CX' = BU    -->   X= -(G^-1)C*X' + (G^-1)*BU
% Let A and R be as shown above, the moments of the MIMO transfer function
% which are the response of the network for impulse excitation are:
% U=[1,1, ...,1]^T --->  X= A*X' + R ---> X= A*sX + R
% To find M0 and M1 then:
% @ DC/ s=0: Mo = X(s=0) ---> X(s=0)=M0 =R;  where R=G^-1* B
%            M1 = X'(s=0)---> G X'(0) = - CX(0) ---> M1 = A * R;  where A=-G^-1* C
%
% and in general:  Mi=A^i*R  or  Mn+1=A*Mn, M0=R
%--------------------------------------------------------------------------
te=tic;
% G*R=B ---> R = (G^-1)B --->M0=R;
R =   Glu.Q* ( Glu.U \( Glu.L  \ (Glu.P * B) ) ) ;
xcput(toc(te), 'FB-to-Form-R', 1, size(Glu.L,1));

end

