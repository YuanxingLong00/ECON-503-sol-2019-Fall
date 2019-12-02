
clear

% Parameters

Par.beta = 0.9;
Par.delta = 0.0196;
Par.rho = 0.9;
Par.sigma = 0.06;
Par.alpha = 0.3;
Par.gamma = 0.5;
Par.K0 = 3;
Par.Z0 = 0;

% Steady State
Zstar = 0;
Kstar = ((1/Par.beta - 1 + Par.delta)./(Par.alpha * exp(Zstar))).^(1/(Par.alpha-1));

% Grid for Capital K
Grid.nK = 200;
Grid.K = linspace(1e-2,2*Kstar, Grid.nK)';        
[Grid.K0d, Grid.k0] = min(abs(Grid.K - Par.K0));    % K0

% Grid for Shocks Z
Par.meanZ = 0;
Grid.nZ = 7;  
[Grid.Z, Grid.PZ]  = mytauchen(Par.meanZ, Par.rho, Par.sigma, Grid.nZ);
[Grid.dZ0, Grid.z0] = min(abs(Grid.Z-Par.Z0));      % Z0
Grid.PZ0 = zeros(1,Grid.nZ);Grid.PZ0(Grid.z0) = 1;

% Product of the two grids
[ZZ,KK] = meshgrid(Grid.Z,Grid.K);
Grid.KK = KK(:);
Grid.ZZ = ZZ(:);



% Initial guess of policy consumption

C = f(Par, Grid.KK,Grid.ZZ ) - Par.delta*Grid.KK;               

b = PolyMatrix(Grid.KK,Grid.ZZ) \ C;

Kp0 = zeros(size(Grid.KK));

MAXIT = 1000;

for it = 1:MAXIT
    
    Kp = f(Par,Grid.KK,Grid.ZZ) - PolyMatrix(Grid.KK,Grid.ZZ)*b;
    
    C = EulerRHS(Par,Grid,Kp,b);
    
    b = PolyMatrix(Grid.KK,Grid.ZZ) \ C;
    
    test = max(abs(Kp - Kp0));
	
    Kp0 = Kp;

    disp(['iteration ' num2str(it) ', test = ' num2str(test)])
    
    if test < 1e-6
       break
    end
    
end

Kp = real(Kp);
bKp = PolyMatrix(Grid.KK,Grid.ZZ) \ Kp;
% KpZ0 = KpMatrix(:,Grid.z0);



% ++++++++++++++++++ Solutions to the Problems ++++++++++++++++

% (3) Plot 3-D Policy Function
KpMatrix = reshape(Kp,Grid.nK,Grid.nZ);
figure1 = surf(KK,ZZ,KpMatrix);
xlabel("Saving");
ylabel("Shock ln(A)");
zlabel("Policy Function");
saveas(figure1,'HW1_2_3');
saveas(figure1,'HW1_2_3.png');

% (4) Plot 3-D Value Function
Vmax = ValueFunc(Par, Grid, Kp);
Vmax = reshape(Vmax,Grid.nK,Grid.nZ);
figure2 = surf(KK,ZZ,Vmax);
xlabel("Saving");
ylabel("Shock ln(A)");
zlabel("Value Function");
saveas(figure2,'HW1_2_4');
saveas(figure2,'HW1_2_4.png');


% (5) Plot Expected Household Capital for next 20 periods
KMatrix = SimMany(Par, Grid, bKp, 20);
EK = mean(KMatrix,2);
t = (0:20)';
figure3 = plot(t,EK);
xlabel("Time");
ylabel("Expected Household Capital");
saveas(figure3,'HW1_2_5');
saveas(figure3,'HW1_2_5.png');

% (6) Find stable distribution and plot it
Tlimit = 400;
KMatrix = SimMany(Par, Grid, bKp, Tlimit);
Klimit = KMatrix(Tlimit,:);
figure4 = histogram(Klimit);
saveas(figure4,'HW1_2_6');
saveas(figure4,'HW1_2_6.png');

% +++++++++++++++ Function Definitions +++++++++++++++++++++++

% ========================
% mytauchen

function [Z, PZ] = mytauchen(e,rho,sigma,N)

m     = 3;
Z     = zeros(N,1);
PZ    = zeros(N,N);
Z(1)  = e - m*sqrt(sigma^2/(1-rho^2));
Z(N)  = e + m*sqrt(sigma^2/(1-rho^2));
step    = (Z(N)-Z(1))/(N-1);

for i=2:(N-1)
   Z(i) = Z(i-1) + step; 
end

for j = 1:N
    for k = 1:N
        if k == 1
            PZ(j,k) = normcdf((Z(1) - (1-rho)*e - rho*Z(j) + step/2) / sigma);
        elseif k == N
            PZ(j,k) = 1 - normcdf((Z(N) - (1-rho)*e - rho*Z(j) - step/2) / sigma);
        else
            PZ(j,k) = normcdf((Z(k) - (1-rho)*e - rho*Z(j) + step/2) / sigma) - ...
                      normcdf((Z(k) - (1-rho)*e - rho*Z(j) - step/2) / sigma);
        end
    end
end

end



% ======= Production Function =======

function f = f(Par, K, Z)

    f = exp(Z) .* K.^Par.alpha + (1-Par.delta)*K;

end

% ======= Derivative of Production Function =======

function d = fprime(Par, K, Z)

    d = exp(Z) .* Par.alpha .* K .^(Par.alpha-1) + (1-Par.delta);

end


% ======================
% Polynomial Matrix [1,A,Y,AY,A^2,Y^2]
% Inputs:
% A: n * 1
% Y: n * 1 

function M = PolyMatrix(K, Z)

M = [ones(size(K)) K Z.*ones(size(K)) K.*(Z.*ones(size(K))) K.^2 (Z.*ones(size(K))).^2];

end


% =================
% RHS of Euler Equation

function C = EulerRHS(Par,Grid,Kp,b)

% inputs
% Kp     nK*Grid.nZ x 1  array of savings
% b      consumption function polynomial coefficients
%
% outputs
% C      consumption implied by the Euler equation


MpOfZ = zeros(Grid.nK*Grid.nZ, Grid.nZ);

for iZp = 1:Grid.nZ
    MpOfZ(:,iZp) = Par.beta * (PolyMatrix(Kp,Grid.Z(iZp)) * b) .^(-Par.gamma) .* fprime(Par,Kp,Grid.Z(iZp));
end

MU = sum( kron(Grid.PZ,ones(Grid.nK,1)) .* MpOfZ ,2);
C = (MU).^(-1/Par.gamma);

end


% =============================
% From policy function to value function

function Vmax = ValueRHS(Par, K, Kp, Z, b)

Vmax = sqrt(f(Par,K,Z) - Kp) + Par.beta * PolyMatrix(Kp, Z) * b;

end

function Vmax = ValueFunc(Par, Grid, Kp)

b = zeros(6,1);
V0 = zeros(size(Grid.KK));

MAXIT = 2000;

for it = 1:MAXIT
    
    % Get Maximized RHS of the Bellman Equation
    Vmax = ValueRHS(Par, Grid.KK, Kp, Grid.ZZ, b);
    
    % Expectated Value of Next Period's Value
    EV = reshape(Vmax, Grid.nK, Grid.nZ) * Grid.PZ';
    
    % Updating p, the linear approximator 
    b = PolyMatrix(Grid.KK, Grid.ZZ) \ EV(:);
    
    test = max(abs(Vmax - V0));
    V0 = Vmax;
    
    disp(['value iteration ' num2str(it) ', test = ' num2str(test)])
    if test < 1e-6
        break
    end
    
end

end

% ==============================
% Simulation

function Sim = SimOnce(Par, Grid, bKp ,T)

% Inputs:
% Par       Parameter structure
% bKp       Polynomial coefficients for polynomial for Kp policy rule


Sim.K = zeros(T+1,1);
Sim.Z = zeros(T+1,1);

Sim.K(1) = Par.K0;
Sim.Z(1) = Par.Z0;
Eps = Par.sigma * randn(T+1,1);


for t = 2:T+1
    Sim.K(t) = PolyMatrix(Sim.K(t-1),Sim.Z(t-1)) * bKp;
    Sim.Z(t) = (1-Par.rho)*Par.meanZ + Par.rho * Sim.Z(t-1) + Eps(t);
    
    [~,k] = min(abs(Grid.K - Sim.K(t)));        % k to keep on the grid
%     disp(["k= " num2str(k)]);
%     disp(["Grid.K= " num2str(Grid.K(k))]);
    Sim.K(t) = Grid.K(k);
    
%     [~,z] = min(abs(Grid.Z - Sim.Z(t)));        % z to keep on the grid
% %     disp(["k= " num2str(z)]);
% %     disp(["Grid.K= " num2str(Grid.Z(z))]);
%     Sim.Z(t) = Grid.Z(z);
end

end

function KMatrix = SimMany(Par, Grid, bKp, T)

MAXSIM = 10000;

KMatrix = zeros(T+1,MAXSIM);

for i_sim = 1:MAXSIM
    
    Sim = SimOnce(Par, Grid, bKp, T);
    
    KMatrix(:, i_sim) = Sim.K;
    
end

end