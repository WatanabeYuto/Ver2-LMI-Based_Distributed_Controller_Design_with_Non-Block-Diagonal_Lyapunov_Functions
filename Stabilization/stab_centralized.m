function [K_opt,P_opt] = stab_centralized(params,flag)

    n = params.n;
    % G = params.G;

    A = params.A;
    B = params.B;

    %%  -------- solve LMI -------- 
    yalmip('clear')

    ops = sdpsettings('solver',params.solver);
    ops.verbose = flag;

    Q = sdpvar(n,n,'symmetric');
    Z = sdpvar(n,n,'full');
    
    % Objective and constraints
    
    LMI = [];
    
    LMI_1 = A*Q + B*Z;
    
    LMI = [LMI_1 + LMI_1' <= 0, Q >= 0];

    optimize(LMI,0,ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z);
    Q_opt = value(Q);
    K_opt = Z_opt*inv(Q_opt);
    P_opt = inv(Q_opt);
    

    fprintf('------------------------\n')
    fprintf('---- Centralized control ---\n')
    % fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    % fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    % fprintf(' max of A+BKs eigval (real part): %8.2e \n', max( real(eig( A + B*K_opt )) ));
    fprintf('----------------------------\n')

end