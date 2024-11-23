function output = check_stab(params,K)
    tmp =  max(real(eig( params.A+params.B*K )));
    % tmp = min((eig(P)));
    if tmp <0
        output =  1;
    else
        output = 0;
    end
end