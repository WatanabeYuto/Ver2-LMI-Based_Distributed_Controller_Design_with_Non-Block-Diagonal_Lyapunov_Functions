function output = check_A_eig(num,A_set)
    % output =  max(real(eig( params.A+params.B*K )));
    output = 0;
    for kk = 1:num
        if max(real(eig(A_set{kk}))) >0
            output = output +1;
        end
    end
end