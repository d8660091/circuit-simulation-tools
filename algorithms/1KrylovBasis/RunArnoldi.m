%--------------------------------------------------------------------------
% aMain_BlkArnoldi.m
% My MIMO Simulator using BLOCK ARNOLDI Algorithm
%--------------------------------------------------------------------------
function [Q,H] = RunArnoldi( G, C, B, Nr_Block_Arnoldi_Moments )
    t0=tic;
    [Glu.L, Glu.U, Glu.P, Glu.Q]= lu( G, 0.1 );
    xcput(toc(t0), 'LU(G)', 1, size(G,1));

    R = FormR(Glu, B);

    [~, Nr.Ports]=size(B);
    % Krylov using Arnoldi Algorithm for (partial) reduction to Hessenberg form
    q = Nr_Block_Arnoldi_Moments * Nr.Ports; % The Order of the Reduced Model or The number of the Column Vector in Q

    Arnoldi_Type = 1;    %Default is 1=BlockArnoldi,  BlockArnoldi can also handel the SISO cases
    t0=tic;
    switch( Arnoldi_Type )
        case 1  %MIMO sys wth reOrthogonalization
            [Q, H] =  BlockArnoldi( Glu, C, R, q);
        case 2  %SISO
            [Q, H] =  Arnoldi( Glu, C, R, q);
    end
    xcput(toc(t0), 'Forming-Q', 1, q);

end %Func
