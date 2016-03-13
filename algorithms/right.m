% Computes right projecting matrix V
% Input : os = original system
% Output: U = left projection matrix
%
function [Q,H] = right(os, q)

    % [Alu.L,Alu.U,Alu.P,Alu.Q]= lu(os.A,0.1);
    % R = Alu.Q*(Alu.U\(Alu.L\(Alu.P * os.B))) ;
    R=os.A\os.B;

    [n, N] = size(R);  %n: MNAsize & N: PortsNr
    NrBlkM = ceil(q/N);  %About to form kr(A, R, NrBlkM)
    q = NrBlkM*N;

    H = zeros(q,q);
    Q = zeros(n,q);

    [Q(:,1:N),~] = qr(full(R) ,0);

    for j = 1: NrBlkM
        %     fprintf('\n Block-Moment # %g',j);
        ci = [(j-1)*N+1 : j*N];  %column of the location the Block in H
        Z = -(Alu.Q*(Alu.U\(Alu.L\((Alu.P*os.E)*Q(:,ci)))));
        for eps = 1:2     %Double Orthogonalization {see the "NOTE1" above}
            %%% Modified Gram-Schmidt orthogonalization:
            for i = 1:j   %Moves along the row in H
                ri = [(j-i)*N+1 : (j-i)*N+N ]; %row of the location the Block in H
                h = Q(:,ri).' * Z;
                Z = Z - Q(:,ri)* h;
                H(ri, ci) = H(ri, ci)+h;
            end
        end
        if i<NrBlkM
            [V,~] = qr(Z,0);
            rii = [1+j*N: (j+1)*N];
            Q(:, rii) = V;
            H(rii, ci)= h;
        end
    end

end %EOF
