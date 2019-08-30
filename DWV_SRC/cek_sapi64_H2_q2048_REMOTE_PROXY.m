%% -----------------------------------------------------------------------
% Solve H2 Problem Using DWAVE
% Created by: suksmono@{stei.itb.ac.id, mdrft.com}
% Date: Sept. 2018
% ------------------------------------------------------------------------
clear all; clc
% Create a local connection
addpath sapi-matlab-client-3.0-win64

%% connect to a server, create a solver
myUrl='https://your.given.url.by_dwave/';
myToken='FILL-IN-YOUR DWAVE TOKEN'
myProxy='http://proxy_usr:proxy_passwd@your.proxy.address:port';
remoteConn = sapiRemoteConnection(myUrl, myToken, myProxy);
%
solver = sapiSolver(remoteConn, 'DW_2000Q_VFYC_2_1');
% Retrieve solver properties from solver
props = sapiSolverProperties(solver);
NQ=props.num_qubits;

%% read hi and Jij from an input file
fname= 'dwave_input\H2_q2048.txt';
h = zeros(NQ, 1);
J=sparse(NQ,NQ);
%w= textscan(fname,'%f');
%% read forst row : NQ N-entries
[nq,NE]= textread(fname,'%f %f',1);
% reast all rows
[vr, vc, vval]= textread(fname,'%f %f %f',NE+1);
for m=1:NE
    row=round(vr(m))+1;
    col=round(vc(m))+1;
    val=vval(m);
%     disp(sprintf('r=%g, c=%g, row=%g, col=%g, val=%f', r,c,row,col, val))
    if (row==col)
        h(row,1)=val;
    else
        J(row,col)=val;
    end
end


%% Solve an ising problem with optional parameters 'num_reads' and 'num_spin_,!reversal_transforms'
aconf = sapiSolveIsing(solver, h, J, 'num_reads', 1000);%, 'num_spin_reversal_transforms', 2)
% display
x=1:length(aconf.energies);
figure(1);
subplot(211);bar(x,aconf.energies,'r')
xlabel('Solution Id'); ylabel('Energy')
% figure(2)
subplot(212);bar(x,aconf.num_occurrences)
xlabel('Solution Id'); ylabel('Occurence')
%% check the solution
[vEn, idE]=min(aconf.energies);
% get the solution vector
sol=aconf.solutions(:,idE);
%% qubits mapping of the problem
% 0 <- q0, q4, q32
% 1 <- q1, q5, q33
% 2 <- q2, q6, q34
% 3 <- q3
disp(sprintf('Minimum energy: %2.3f', vEn));
disp(sprintf('H2 and D2'));
vH2=[sol(1), sol(2), sol(3), sol(4) ];
H2=reshape(vH2,2,2);
D2=H2*transpose(H2);
H2,
D2,

