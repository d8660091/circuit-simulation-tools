% Use eigen decomposition for solving lyapunov equation
function P=eig_dec_lyap(A,E)
    A_bar=inv(A)*E;
    [V_bar,D_bar]=eig(full(A_bar));
    P_bar=inv(V_bar')*inv(V_bar);
    P=inv(A')*P_bar*inv(A);
end
