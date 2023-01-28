% Extract solutions from the moment matrix if the flatness condition is
% satisfied.

% Input
% Mmatrix: moment matrix
% x: Yalmip variables
% v: monomial basis

% Output
% sol: a matrix consisting of optimal solutions

function sol = extract_solution(Mmatrix, x, v)
    n = length(x);
    [U, pivot] = rref(Mmatrix', 1e-3);
    rk = length(pivot);
    U = U(1:rk,:)';
    w = v(pivot);
    for i = 1:n
        xw = x(i)*w;
        k = [];
        for j = 1:rk
            k = [k; find(ismember(xw(j),v))];           
        end
        N{i} = U(k,:);
    end
    
    % Create random convex combination
    rands = rand(n,1);
    rands = rands/sum(rands);
    M = zeros(rk);
    for i = 1:n
        M = M + rands(i)*N{i};
    end

    [L,~] = schur(M);
    sol = zeros(n,rk);
    for i = 1:rk
        for j = 1:n
            sol(j,i) = L(:,i)'*N{j}*L(:,i);
        end
    end
end