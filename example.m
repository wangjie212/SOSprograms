spotpath = '../../../Programs/spotless';
addpath(genpath(spotpath));
sdptpath = '../../../Programs/SDPT3-4.0';
addpath(genpath(sdptpath));

% x = msspoly('x', 2);
% a = [1;1];
% f = [x(1)+x(2); x(1) - x(1)^2; x(2) - x(2)^2];
% g = [x(1)*x(2)];
% h = [x(1)^2+x(2)^2-1];
% [blk, At, C, b] = sosprogram(a, f, g, h, x, 2);
% sqlp(blk, At, C, b);

% Extract solutions
% yalmippath = '../../../Programs/YALMIP';
% addpath(genpath(yalmippath));

% sdpvar x1 x2
% obj = x1^4-x1^2+x2^2-2*x1*x2+1;
% [info,sol,data] = solvemoment([],obj,[],3);
% 
% nsol = extract_solution(data.moment{3}, data.x, data.monomials);

% An example with sparsity
x = msspoly('x', 6);
a = [-1];
f = [1+x(1)^4+x(2)^4+x(3)^4+x(4)^4+x(5)^4+x(6)^4+x(1)*x(2)*x(3)+x(3)*x(4)*x(5)+x(3)*x(4)*x(6)+x(3)*x(5)*x(6)+x(4)*x(5)*x(6); -1];
g = [1-x(1)^2-x(2)^2-x(3)^2];
h = [1-x(3)^2-x(4)^2-x(5)^2-x(6)^2];
cliques{1} = [1;2;3];
cliques{2} = [3;4;5;6];
I{1} = [1];
I{2} = [];
J{1} = [];
J{2} = [1];
order = [2;2];
[blk, At, C, b] = sparsesos(a, f, g, h, x, cliques, I, J, order);
sqlp(blk, At, C, b);
