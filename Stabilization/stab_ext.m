function [K_opt,P_opt] = stab_ext(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;

    %%  -------- solve LMI -------- 
    yalmip('clear')

    ops = sdpsettings('solver',params.solver);
    ops.verbose = flag;

    Z = [];
    Q = sdpvar(n,n,'symmetric');
    Slk = diag(sdpvar(n,1));
    % gamma = 1.0981;
    
    L = laplacian(G);
    for i = 1:n
        tmp = [];
        for j = 1:n
            if L(i,j) ~= 0
                tmp = [tmp sdpvar(1,1)];
            else
                tmp = [tmp 0];
            end
        end
        Z = [Z;tmp];
    end
    
    % Objective and constraints
     
    LMI = [];
    
    LMI_1 = kron([0, 1; 0, 0],Q) + [Slk'*A' + Z'*B', Slk'*A' + Z'*B'; -Slk', -Slk'];
    
    LMI = [LMI_1 + LMI_1' <= 0, Q >= 0];


    optimize(LMI,0,ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = fillmissing(value(Z),"constant",0);
    Q_opt = value(Q);
    Slk_opt = value(Slk);
    K_opt = Z_opt*inv(Slk_opt);
    P_opt = inv(Slk_opt')*Q_opt*inv(Slk_opt);
    % eig(value(LMI_1))

    fprintf('------------------------\n')
    fprintf('----Extended LMI---\n')
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    % fprintf(' Norm of K                      : %8.2e \n', norm(K_opt)); 
    fprintf(' max of A+BKs eigval (real part): %8.2e \n', max( real(eig( A + B*K_opt )) ));
    fprintf('--------------------------------\n')

end