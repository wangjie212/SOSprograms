clc
clear all
close all
% restoredefaultpath

% mosekpath   = "../../../mosek";
% sdpt3path   = "../SDPT3";
% sostools    = "../SOSTOOLS";
% spotpath    = "../CBF_framework/matlab/spotless";
% addpath(genpath(mosekpath))

kappa = 3;
umax  = 1;
thmax = 2;

%% A SOS program using SOSTOOLS, MOSEK solves to global optimality
% addpath(genpath(sostools))
x   = mpvar('x',[2,1]);
x1  = x(1);
x2  = x(2);
z   = mpvar('z',[1,1]);
th  = mpvar('th',[1,1]);
var = [x;z;th];

% solve the sos program
% max sum_beta lambda_beta * gamma_beta
% s.t., phi - (sum_beta lambda_beta * theta^beta) = ...
% sigma0*s0 + sigma1*s1 + ... + mu1*h1 + mu2*h2 + ...
% sigma0, sigma1 ... are SOS
% mu1, mu2, ... are polynomials

phi = -1 * x2^2 * (1-x2^2) + 2*umax*z;

% equality constraints
h   = [
    z^2 - x1^2*x2^2;
    x1^2 + x2^2 - th
    ];
% inequality constraints
s   = [
    z;
    th*(thmax - th)
    ];

dmax = max([phi.maxdeg,s.maxdeg,h.maxdeg]);
fprintf("dmax = %d, kappa = %d.\n",dmax,kappa);
assert(2*kappa >= dmax,'increase relaxation order')

% theta monomials
th_basis = monomials(th,0:(2*kappa));
p = th_basis.coefficient * th_basis.degmat;
gamma = -1 * (1 ./ (p+1)) .* (thmax .^ (p));

[lamstar,out] = mysosprogram(gamma, [phi;th_basis], s, h, var, kappa);

%% Same SOS program using your code, MOSEK/SDPT3 returns infeasible
% rmpath(genpath(sostools))
% addpath(genpath(pwd))
% addpath(genpath(spotpath))
% addpath(genpath(sdpt3path))

x   = msspoly('x',2);
x1  = x(1);
x2  = x(2);
z   = msspoly('z',1);
th  = msspoly('th',1);
var = [x;z;th];

phi = -1 * x2^2 * (1-x2^2) + 2*umax*z;

% equality constraints
h   = [
    z^2 - x1^2*x2^2;
    x1^2 + x2^2 - th
    ];
% inequality constraints
s   = [
    z;
    th*(thmax - th)
    ];

% theta monomials
th_basis = monomials(th,0:(2*kappa));
p = 0:(2*kappa);
gamma = -1 * (1 ./ (p+1)) .* (thmax .^ (p));

[blk, At, C, b] = sosprogram(gamma, [phi;-1*th_basis], s, h, var, kappa);
sqlp(blk, At, C, b);





