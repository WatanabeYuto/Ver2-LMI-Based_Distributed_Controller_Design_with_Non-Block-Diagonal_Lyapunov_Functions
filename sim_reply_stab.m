clear all;

%% define parameters
%%%%%%%%% --------------- start ---------------
parameters

A_set = load("A_matrices_stab.mat");
A_set = A_set.A_set;
B_set = load("B_matrices_stab.mat");
B_set = B_set.B_set; 

Graph_set = {};
Graph_set{1} = generate_ringgraph(params.n);
Graph_set{2} = generate_wheelgraph(params.n);

% plot(params.G); 

stab_result = {};
stab_result{1} = [];
stab_result{2} = [];

%%%%%%%%% --------------  end  ---------------

params.solver = 'sdpt3';
% params.solver = 'mosek';
% params.solver = 'sedumi';

for gr = 1:2
    
    params.G = Graph_set{gr};

    for kk = 1:number_of_data
        params.A = A_set{kk};
        params.B = B_set{kk};

        params.C = eye(params.n);
        params.D = eye(params.n);
        %% exogenous noises
        params.Bw = eye(params.n);
        params.Dw = eye(params.n);
        
        %% run 
        %%%%%%%%% --------------- start ---------------

        % proposed method 1
        [K_opt_proposed1,P_opt_proposed1] = stab_proposed1(params,0);

        % proposed method 3
        [K_opt_proposed3,P_opt_proposed3] = stab_proposed3(params,0);

        %combined method
        [K_opt_c,P_opt_c] = stab_combined(params,0);

        %proposed method 2
        [K_opt_2,P_opt_2] = stab_proposed2(params,0);

        % extended LMI
        [K_opt_ext,P_opt_ext] = stab_ext(params,0);

        % SI
        [K_opt_SI,P_opt_SI] = stab_SI(params,0);

        % block-diagonal relaxation:
        [K_opt_diag,P_opt_diag] = stab_diag(params,0);

        % check_stab(params,K_opt_proposed_0_v2)
        stab_chek_kk = [check_stab(params,K_opt_proposed1),check_stab(params,K_opt_proposed3),check_stab(params,K_opt_c),check_stab(params,K_opt_2),check_stab(params,K_opt_ext),check_stab(params,K_opt_SI),check_stab(params,K_opt_diag)];
        stab_result{gr} = [stab_result{gr}; stab_chek_kk];

    end
end

function output = check_stab(params,K)
    tmp =  max(real(eig( params.A+params.B*K )));
    % tmp = min((eig(P)));
    if tmp <0
        output =  1;
    else
        output = 0;
    end
end
