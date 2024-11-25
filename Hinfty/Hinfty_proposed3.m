function [gamma_opt,K_opt,P_opt] = Hinfty_proposed3(params,flag)

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
    N = E*inv(E'*E)*E'; 

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

    Z_e = [];
    Q_e = [];
    gamma = sdpvar(1,1);
    
    for l = 1:length(cliques)
        tmp_q = sdpvar(size(EE{l},1),size(EE{l},1),'symmetric');
        tmp_z = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        Q_e = blkdiag(Q_e, tmp_q);
        Z_e = blkdiag(Z_e, tmp_z);
    end

    dim = size(Q_e,1);

    LMI = [];

    LMI_s = A_e * Q_e + B_e * Z_e;
    
    LMI_1 = [LMI_s + LMI_s', Bw_e, Q_e*C_e'+ Z_e'*D_e';
            Bw_e', - gamma * eye(n), Dw';
            C_e*Q_e + D_e*Z_e, Dw, -gamma * eye(n)];
    
    rho = sdpvar(1,1);
 
    Theta = blkdiag(rho*M,zeros(n),zeros(n));
    LMI = [LMI_1 + Theta <= -1/10^8*eye(dim+n+n), Q_e >=1/10^8*eye(dim)];
    
    optimize(LMI,gamma,ops)
    
    
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z_e);
    Q_opt = value(Q_e);
    K_opt = inv(E'*E)*E'*(Z_opt*inv(Q_opt))*E;
    gamma_opt =  value(gamma);
    P_opt = E' * inv(Q_opt) * E;


    fprintf('------------------------\n')
    fprintf('----Proposed method 3---\n')
    fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (real part): %8.2e \n', max( real(eig( A + B*K_opt )) ));
    fprintf(' max of PAt+APs eigval (real part): %8.2e \n', max( real(eig( P_opt*(A + B*K_opt) +  (P_opt*(A + B*K_opt))' )) ));
    fprintf('------------------------\n')

end