parameters;

number_of_data = 100;

A_set = {};
B_set = {};


max_it = 10;
for j = 1:number_of_data
    
    for jj = 1:max_it 
        A_set{j} = randn(params.n,params.n);
        B_set{j} = diag([ones(params.n-dd,1);zeros(dd,1)]);

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

