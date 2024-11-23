function [gamma_opt,K_opt,P_opt] = Hinfty_centralized(params,flag)

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

    % %% dilation
    % A_e = dilation_by_E(A,E);
    % B_e = dilation_by_E(B,E);
    % C_e = C * inv(E'*E) * E';
    % D_e = D * inv(E'*E) * E';
    % Bw_e = E * Bw;
    % % Dw_e = Dw;
    
    % cliques = maximalCliques(adjacency(G)); 

    %%  -------- solve LMI -------- 
    yalmip('clear')

    ops = sdpsettings('solver',params.solver);
    ops.verbose = flag;

    Q = sdpvar(n,n,'symmetric');
    Z = sdpvar(n,n,'full');
    gamma= sdpvar(1,1);
    
    % Objective and constraints
    
    LMI = [];
    
    LMI_1 = A*Q + B*Z;
    
    LMI_1 = [LMI_1 + LMI_1', Bw, Q*C'+ Z'*D';
            Bw', - gamma * eye(n), Dw';
            C*Q + D*Z, Dw, -gamma * eye(n)];
    
    LMI = [LMI_1 <= 0, Q >= 1/10^8*inv(E'*E)];

    optimize(LMI, gamma,ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z);
    Q_opt = value(Q);
    K_opt = Z_opt*inv(Q_opt);
    gamma_opt=value(gamma);
    P_opt = inv(Q_opt);
    

    fprintf('------------------------\n')
    fprintf('---- Centralized control ---\n')
    fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (real part): %8.2e \n', max( real(eig( A + B*K_opt )) ));
    fprintf('----------------------------\n')

end