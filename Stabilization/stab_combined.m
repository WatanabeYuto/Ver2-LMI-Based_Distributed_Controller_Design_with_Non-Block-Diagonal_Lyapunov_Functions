function [K_opt,P_opt] = stab_combined(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;

    Bw = params.Bw;
    Dw = params.Dw;

    E = generate_Ematrix(n,G);
    EE = generate_Ematrix_cell(n,G);
    M = eye(size(E,1)) - E*inv(E'*E)*E';

    %% dilation
    A_e = dilation_by_E(A,E);
    B_e = dilation_by_E(B,E);
    C_e = C * inv(E'*E) * E';
    D_e = D * inv(E'*E) * E';
    Bw_e = E * Bw;
    % Dw_e = Dw;
    
    cliques = maximalCliques(adjacency(G)); 

    %%  -------- solve LMI -------- 
    yalmip('clear')

    ops = sdpsettings('solver',params.solver);
    ops.verbose = flag;

    Q = sdpvar(n,n,'symmetric'); 
    mu = sdpvar(1,1);  
    Q_e = mu * eye(size(E,1)) + E*Q*E';
    Z_e = [];
    G_e = [];
    
    for l = 1:length(cliques)
        tmp_q = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        tmp_z = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        G_e = blkdiag(G_e, tmp_q);
        Z_e = blkdiag(Z_e, tmp_z);
    end
    
    dim = size(Q_e,1);

    LMI = [];

    LMI_1 = kron([0, 1; 0, 0],Q_e) + [G_e'*A_e' + Z_e'*B_e', G_e'*A_e' + Z_e'*B_e'; - G_e', -G_e'];
    
    rho = sdpvar(1,1);
    eta = sdpvar(1,1); 
   
    LMI = [LMI_1 + LMI_1' + blkdiag(rho*M,rho*M) <= 0, Q_e >=0, G_e'*M + M*G_e >= eta*M, eta >=0];


    
    optimize(LMI,0,ops)
    
    %%  --------  end LMI  --------


    %% results 
    Z_opt = fillmissing(value(Z_e),"constant",0);
    Q_opt = value(Q_e);
    G_e_opt =value(G_e);
    K_opt = inv(E'*E)*E'*(Z_opt*inv(G_e_opt))*E;
    % P_opt = inv(E' * Q_opt * E);
    P_opt = E' * inv(G_e_opt')*Q_opt* inv(G_e_opt) * E;
    value(rho/eta)



    fprintf('------------------------\n')
    fprintf('----combined method---\n')
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    % fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (real part): %8.2e \n', max( real(eig( A + B*K_opt )) ));
    fprintf('------------------------\n')

end