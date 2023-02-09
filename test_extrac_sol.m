clc
clear
close all
restoredefaultpath

mosekpath   = "../../../mosek";
spotpath    = "../CBF_framework/matlab/spotless";
addpath(genpath(mosekpath))
addpath(genpath(spotpath))
addpath(genpath(pwd))

load("univariate.mat")

kappa = deg(vnu)/2 + 1;

thmax = 2.0;
problem.vars            = th; % univariate 
problem.objective       = -1 * vnu; % goal is to maximize vnu
problem.inequality      = [th*(thmax-th)]; % theta is between 0 and 2
[SDP,info]              = dense_sdp_relax(problem, kappa);
prob = convert_sedumi2mosek(SDP.sedumi.At,...
                            SDP.sedumi.b,...
                            SDP.sedumi.c,...
                            SDP.sedumi.K);
[~,res]              = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
Xmom = Xopt{1};
disp(rank(Xmom,1e-3))
disp(rank(Xmom(1:end-1,1:end-1),1e-3))

sol = extract_solution(Xopt{1},th,info.v);

disp(sol)

