%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---- Finding Q:---------------------------------------------------------
%       Determining the orthogonal basis Q of the Krylov subspace
%       Kr(A,R,n) where n<m.
%--------------------------------------------------------------------------
%%
function [Q, H] = Arnoldi(G, C, R, q)
    R = full(R);
    Q = R/norm(R,2);
    for nc = 1: q %q is No_of_column_in_Q
        Z = -G\(C*Q(:,nc));
        for nr = 1:nc %nr=Q_row_Count
            H(nr,nc) = Q(:,nr).'* Z;
            Z = Z - H(nr,nc)*Q(:,nr);
        end

        if nc < q
            H(nc+1,nc) = norm(Z,2);
            if H(nc+1,nc)==0; break; end
            Q(:,nc+1 ) = Z / norm(Z,2);
        end
    end

end %EOF
