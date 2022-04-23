clear; format short e

%% State Space Matricies
%dymaics matrix
A = [1.007  0.0050;
    3.0977 1.0770];
%input matrix
B= [-0.0186;
    -7.4626];
%output matrix
C = [0 1];
%feedthrough matrix
D = 0;

pts = 100;
time = linspace(.01,2,pts);

%% LQR penelty matricies
plant_lqr = ss(A,B,C,D);
Q = [1000 0;
    0 1];
R = 1;
P = [1 0;   %nominal terminal cost
    0 1];

%% Terminal cost for stability
Alen = length(A);
Blen = size(B);
Blen = Blen(2);

%variables to solve for
P = sdpvar(Alen,Alen); %terminal cost matrix
K = sdpvar(Blen,Alen); %state feedback control matrix

%Symmetric LMI found using schur compliment
LMI = [ -P              P*A'+K'*B'          P                   K';
        (P*A'+K'*B')'   -P                  zeros(Alen,Alen)    zeros(Alen,Blen);
        P'              zeros(Alen,Alen)    -inv(Q)             zeros(Alen,Blen)
        K               zeros(Blen,Alen)    zeros(Blen,Alen)    -inv(R)];
    
const = [LMI<=0, P>=0]; %constraints
solvesdp(const); %solver

%Convert solutions to doubles
P=double(P); %P_tilde
K=double(K); %K_tilde

%Solve for P and K in terms of P_tilde and K_tilde
P = inv(P); %P = inv(P_tilde)
K = K*P; %K = K_tilde*P

% P = [9.5146e+03   4.8747e+01;   %terminal cost for stability  
%      4.8747e+01   2.2288e+00];

%% Batch Matricies
N = [1,5,10]; % horizon lengths
for Ns = 1:length(N)
    Np = num2str(N(Ns));
    % create batch matricies for each horizon length
    % batch matrix program written by Leila Bridgeman
    eval(strcat('[Ab_',Np,',Bb_',Np,',Qb_',Np,',Rb_',Np,...
        ']=URHC_QO_BatchMatrices(A,B,Q,R,P,',Np,')'));
end

%% Postive invariant set for terminal region; ensures feasability
system = LTISystem('A', A, 'B', B);
system.x.min = [-3; -50]; %lower state bounds (pos, vel)
system.x.max = [3; 50]; %upper state bounds
system.u.min = -10; %input lower bound
system.u.max = 10; %input upper bound
InvSet = system.invariantSet(); %find invariant set
InvSet.computeHRep; %write in H representation [H_x | K_x]
Hrepresent = InvSet.H; 
Hxf = Hrepresent(:,1:2);
Kxf = Hrepresent(:,3);


%% Implementing Controller
% zero matrices: rows = hoizon length, columns = time step
pos = zeros(3,pts); % position x
vel = zeros(3,pts); % velocity xdot
input = zeros(3,pts); % input u

%cycle through random initial conditions
for k=1:1:1
    
    % random initial conditions within range of constraints
    rng(390); %random seed
    x0 = [randi([-3,3]); randi([-50,50])];
    x0 = [3,30]';
    
    tic % start stopwatch

    %cycle through horizon lengths
    for Ns = 1:length(N)
        x=x0; 
        Np = num2str(N(Ns));
        
        %build constraint matricies
        Hx = zeros(4*N(Ns)+6,2*(N(Ns)+1));
        Kx = zeros(4*N(Ns)+6,1);
        Hu = zeros(2*N(Ns),N(Ns));
        Ku = 10*ones(2*N(Ns),1);
        for z = 0:(N(Ns)-1)
            Hx(1+4*z,1+2*z) = 1;
            Hx(2+4*z,1+2*z) = -1;
            Hx(3+4*z,2+2*z) = 1;
            Hx(4+4*z,2+2*z) = -1;
            Kx(1+4*z) = 3;
            Kx(2+4*z) = 3;
            Kx(3+4*z) = 50;
            Kx(4+4*z) = 50;
            Hu(1+2*z,1+z) = 1;
            Hu(2+2*z,1+z) = -1;
        end
        Hxsz = size(Hx);
        Hfsz = size(Hxf);
        Hx(Hxsz(1) - Hfsz(1)+1:Hxsz(1), Hxsz(2) - Hfsz(2)+1:Hxsz(2)) = Hxf;
        Kx(length(Kx) - length(Kxf)+1:length(Kx)) = Kxf;
        
        %cycle through time steps
        for i=1:length(time)
            yalmip('clear')
            
            %u batch
            u_bar = sdpvar(N(Ns),1);
            
            %x batch
            eval(strcat('x_bar = Ab_',Np,'*x+Bb_',Np,'*u_bar;'));
            
            %Objective function x'Qx+u'Ru
            eval(strcat('objective = x_bar''','*Qb_',Np,...
                '*x_bar + u_bar''','*Rb_',Np,'*u_bar;'));

            constraint = [Hx*x_bar<=Kx, Hu*u_bar<=Ku];
            constraint = [];
        
            %optimize objective function
            sol=optimize(constraint,objective);
            
            %Store only the first input %%%%%%%%%%%%try graphing all
            %trajectories for fun
            u = value(u_bar(1));
            
            %Storing pos vel and inputs
            pos(Ns,i) = x(1);
            vel(Ns,i) = x(2);
            input(Ns,i) = u;
            
            %Calculating the next x with optimized input
            x = A*x+B*u;
        end
    end
    toc
end

%% Create contraint vectors for plotting
pos_con_n = ones(1,100)*3;
vel_con_n = ones(1,100)*50;
vol_con_n = ones(1,100)*10;

%% Position Plot
figure()
hold on
plot(time,pos)
legend({'N=1','N=5','N=10'})
title('Ball Position Vs. Time')
ylabel('Position(mm)')
xlabel('Time(s)')
patch([time fliplr(time)], [pos_con_n min([-3 3])*ones(size(pos_con_n))], 'g')  
alpha(.09)
legend({'N=1','N=5','N=10','pos constraint'})
hold off
legend({'N=1','N=5','N=10','Viable Region'})
legend({'N=1','N=5','N=10','Constrained Region'})
ylim([-4,4])

%% Velocity Plot
figure()
hold on
plot(time,vel)
legend({'N=1','N=5','N=10'})
title('Ball Velocity Vs. Time')
ylabel('Velocity (mm/s)')
xlabel('Time(s)')
patch([time fliplr(time)], [vel_con_n min([-50 50])*ones(size(vel_con_n))], 'g')  
alpha(.09)
hold off
legend({'N=1','N=5','N=10','Viable Region'})
legend({'N=1','N=5','N=10','Constrained Region'})
ylim([-100,60])

%% Input plot
figure()
hold on
plot(time,input)
legend({'N=1','N=5','N=10'})
title('Maglev Voltage Vs. Time')
patch([time fliplr(time)], [vol_con_n min([-10 10])*ones(size(vol_con_n))], 'g')  
ylabel('Input Voltage (V)')
xlabel('Time(s)')
alpha(.09)
hold off
legend({'N=1','N=5','N=10','Viable Region'})
legend({'N=1','N=5','N=10','Constrained Region'})
ylim([-11,11])

%% Batch Matrix Generator
function [Ab,Bb,Qb,Rb] = URHC_QO_BatchMatrices(A,B,Q,R,P,N)
% function [Ab,Bb,Qb,Rb] = URHC_QO_BatchMatrices(A,B,Q,R,P,N)
%
% This function computes the batch cost and dynamics matrices for
% unconstrained receeding horizon control with a quadratic objective and
% LTI dynamics.
%
% Inputs
%   A,B   = discrete time dynamics and input matrix for plant
%   Q,R   = running cost, x'Qx + u'Ru
%   P     = terminal cost x'Px
%   N     = horizon length
%
% Outputs
%   Ab,Bb = batch dynamics matrices designed so that
%           xb = Ab x(t) + Bb ub, where
%               xb = [x_0'; x_1';...;x_N'],
%               ub = [u_0'; u_1';...;u_(N-1)'],
%               x_0 = x(t), and
%               x_(i+1) = A x_i + B u_i.
%  Qb,Rb = batch cost matrices designed so that
%           xb'Qbxb + ub'Rub = x_N'Px_N + sum_(i=0)^(N-1) ( x_i'Qx_i + u_i'Ru_i )

n = size(A,1); % # states
m = size(B,2); % # inputs

Qb = blkdiag(kron(eye(N),Q),P);
Rb = kron(eye(N),R);

Ab = zeros((N+1)*n,n); % pre-allocate space for Ab
Ab(1:n,:) = eye(n); % initialization
Bb = zeros((N+1)*n,N*m); % pre-allocate space for Bb

for iter = 1:N
    ri = iter*n + 1; % first row in Ab and Bb to be edited
    rf = (iter+1)*n; % last row in Ab to be edited
    Ab(ri:rf,:) = A*Ab(ri-n:ri-1,:); % add entry to a 
    Bb(ri:end,1:end-(iter-1)*m) = Bb(ri:end,1:end-(iter-1)*m)...
                                + kron(eye(N-iter+1),Ab(ri-n:ri-1,:)*B);
end
end