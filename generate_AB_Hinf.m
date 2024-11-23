
A_set = {};
B_set = {};


max_it = 5000; % iterate sufficiently many timez to obtain a stabilizable plant
dd = 1;
for j = 1:50
    
    for jj = 1:max_it 
        eig_A_real = - rand(params.n/2,1)/2;

        tmp = randperm(params.n/2);
        eig_A_real(tmp(1:3)) = - eig_A_real(tmp(1:3)); % there are three unstable eigenvalues


        eig_A_real = [eig_A_real; eig_A_real];
        eig_A_complex = rand(params.n/2,1);
        eig_A_complex = [eig_A_complex; - eig_A_complex];
        eig_A = complex(eig_A_real, eig_A_complex);

        tmp = rand(params.n,params.n);
        AA = tmp - place(tmp,eye(params.n), eig_A); 
        tmp = rand(params.n,params.n);
        A_set{j} = tmp \ AA * tmp;

        
        B_set{j} = eye(params.n);
        B_set{j}(1,1) = 0;
        

        %% stabilizability
        eig_A = eig(A_set{j});
        tmp = 0;
        for ii = 1:params.n
            if real(eig_A(ii)) >0 && rank([B_set{j}, A_set{j}-eig_A(ii)*eye(params.n)]) < params.n
                tmp = 1;
                break
            end
        end

        if tmp == 0
            break
        end
    end

end

save 'A_matrices_Hinf.mat' A_set
save 'B_matrices_Hinf.mat' B_set