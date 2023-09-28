function [ustar,out] = mysosprogram(a, f, g, h, x, d)

% make sure column vectors
a = a(:);
f = f(:);
g = g(:);
h = h(:);
x = x(:);

% initialize SOS program
% dim_u   = length(f) - 1;
% u       = dpvar(mpvar('u',[dim_u,1]));

prog        = sosprogram(x);

u_name = {};
for i = 1:length(f)-1
    u_name{i} = sprintf('u_%d',i);
end
u = dpvar(u_name);
u = u(:);
prog = sosdecvar(prog,u);

v = u'*f(2:end);

% [prog,v]    = sospolyvar(prog,f(2:end));

% multipliers for equality constraints
t = [];
for i = 1:length(h)
    hi      = h(i);
    di      = hi.maxdeg;
%     di      = polynomialDegree(hi);
    ti_basis    = monomials(x,0:(2*d-di));
    [prog,ti]   = sospolyvar(prog,ti_basis);
    t = [t;ti];
end

% multipliers for inequality constraints
g = [monomials(x,0);g];
s = [];
for i = 1:length(g)
    gi      = g(i);
    di      = gi.maxdeg;
%     di      = polynomialDegree(gi);
    sigi    = floor((2*d - di)/2);
    si_basis    = monomials(x,0:sigi);
    [prog,si]   = sossosvar(prog,si_basis);
    s           = [s;si];
end

eq   = f(1) - v - t'*h - s'*g;
prog = soseq(prog,eq);
prog = sossetobj(prog,a'*u);

options.solver = 'mosek';
prog = sossolve(prog,options);

ustar = double(sosgetsol(prog,u));
out   = prog;
end