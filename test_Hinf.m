%% run 

gamma_result = [];
stab_result_Hinf = [];


%%%%%%%%% --------------- start ---------------

% % H infty
% proposed method:
[gamma_opt_proposed1,K_opt_proposed1,P_opt_proposed1] = Hinfty_proposed1(params,0);

% proposed method v2:
[gamma_opt_proposed2,K_opt_proposed2,P_opt_proposed2] = Hinfty_proposed3(params,0);

% proposed method with rho=:
[gamma_opt_proposed3,K_opt_proposed3,P_opt_proposed3] = Hinfty_proposed2(params,0);

% SI
[gamma_opt_SI,K_opt_SI,P_opt_SI] = Hinfty_SI(params,0);

% block-diagonal relaxation:
[gamma_opt_diag,K_opt_diag,P_opt_diag] = Hinfty_diag(params,0);

% centralized controllr:
[gamma_opt_cen,K_opt_cen,P_opt_cen] = Hinfty_centralized(params,0);


gamma_result = [gamma_opt_proposed1,gamma_opt_proposed2,gamma_opt_proposed3,gamma_opt_SI,gamma_opt_diag,gamma_opt_cen];
stab_result_Hinf = [check_stab_A(params,K_opt_proposed1),check_stab_A(params,K_opt_proposed2),check_stab_A(params,K_opt_proposed3),check_stab_A(params,K_opt_SI),check_stab_A(params,K_opt_diag),check_stab_A(params,K_opt_cen)];

function output = check_stab_A(params,K)
    output =  max(real(eig( params.A+params.B*K )));
end