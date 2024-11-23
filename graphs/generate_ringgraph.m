function G = generate_ringgraph(n)

    % arguments
    %     prob = 0.3 % probability to generate edges
    % end
    
    Adj = zeros(n,n); % adjacency matrix
    
    for i = 1:n-1
        Adj(i,i+1) = 1;
        Adj(i+1,i) = 1;
    end

    Adj(n,1) = 1;
    Adj(1,n) = 1;

    
    G = graph(Adj);
end