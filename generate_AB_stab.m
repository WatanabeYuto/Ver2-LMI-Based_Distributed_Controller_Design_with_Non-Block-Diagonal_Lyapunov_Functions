
A_set = {};
B_set = {};

max_it = 5000;
dd = 1;
for j = 1:200
    
    for jj = 1:max_it 
        A_set{j} = randn(params.n,params.n);

        B_set{j} = eye(params.n);
        B_set{j}(1,1) = 0;
        B_set{j}(round(params.n/2),round(params.n/2)) = 0;
        

        % check the stabilizability
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

save 'A_matrices_stab.mat' A_set
save 'B_matrices_stab.mat' B_set