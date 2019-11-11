%% Question 2 
%% (1) 
% The sequenc problem is 
%
% $$ \max_{ c_t,k_{t+1}} \sum_{t=0}^ \infty \beta U(c_t) $$ subject to 
% $$ c_t + k_{t+1} \leq A_t k_t^\alpha  +(1-\delta)k_t $$ ; $$ ln A_{t+1} =\rho ln A_{t} +u_t $$; $$ c_t \geq 0
% $$; $$ k_{t+1} \geq 0 $$. 
%
% The recursive formulation is 
%
% $$ V(k,A) = \max_{ k'\in [0, Ak^\alpha +(1-\delta)k ],} U(c)
% + \beta E \big( V(k',A')| A \big) $$. 
%
% State variables are $k_t, A_t$. Control variables are $k_{t+1}, c_t$.  
%% (2) 
clear
gamma=0.5;
alpha=0.3;
beta=0.9;
delta=0.95;
a0=3;
A0=1;
rho=0.9;
NA= 9;
[prob,grid_lnA,invdist]=tauchen(NA,0,rho,1);
grid_A= exp(grid_lnA); %% grid set for productivity
xe= (1-beta+delta*beta)/(alpha*beta);
ke= xe^(1/(alpha-1)); %% Steady state capital when A=1 is a constant. 
grid_k = 0.02:0.02:5 ; %% grid set for capital
Nk=size(grid_k,2);
gc= ones(Nk,NA).*0.2; % policy function for consumption
gc_new=ones(Nk,NA).*0.2;
gc1=ones(Nk,NA).*0.2;
gk_ind=ones(Nk,NA); % index policy function for capital
gk=ones(Nk,NA);
err= 10^(-4);
iter=1000;
for i=1:iter
    for knew=1:Nk
        for a=1:NA
            X0= [0.5,0.25];
            X=runnested(knew, a, prob, gc, grid_A,grid_k, X0);
            gc1(knew,a)= X(1); 
            k=int64( find(grid_k > X(2),1,'first') );
            if k ~= 0
                k = Nk;
            end
            gk_ind(k,a)=knew;
        end        
    end
    gk_ind =int64 (gk_ind); % gk_new is index.
    gk = grid_k(gk_ind);
    for k=1:Nk
        for a=1:NA
            gc_new(k,a)= min( max( gc1( gk_ind(k,a), a) ,0), grid_A(a)* grid_k(k)^alpha +(1-delta)*grid_k(k));
        end
    end
    if norm(gc-gc_new)< err 
           fprintf('Number of iteratios: %d, Convergence Achieved. Used tolerance= %1.6f . \n', i, err);
         break
    else
         gc=gc_new;
    end
end
save('EGM','grid_A','grid_k','gc','gk','gk_ind','prob','invdist','Na','Nk');

%% (3)
clear 
load('EGM','grid_A','grid_k','gc','gk','gk_ind','prob','invdist','Na','Nk');
gridk=repmat(grid_k',1,Na);
gridA=repmat(grid_A,Nk,1);
mesh(gridk,gridA,gk)

%% (4) 
% Observe that the consumption is bounded when state space is bounded.
% Therefore, the utility is also bounded. We can approximate the value
% function by computing previous N period total utility. 
N=1001;
B=1000;
V=zeros(Nk,NA,B);
 for b=1:B
    for k=1:Nk
         for a=1:NA
             A=zeros(N,1);
             A(1)=a;
             K=zeros(N,1);
             K(1)=k;
             for j=2:N
                 A(j)= Gen_yind(prob,A(j-1));
                 K(j)= gk_ind(K(j-1),A(j-1));
             end
             
             s=0;
             for j=1:N-1
                 s = s+ beta^j *( grid_A(A(j))*grid_k(K(j))^alpha + (1-delta)* grid_k(K(j))-grid_k(K(j+1)) )^gamma;
             end
             V(k,a,b) = s;
         end
    end
 end
 meanV=mean(V,3);
 ValueFun= reshape(meanV,Nk,NA);
 mesh(gridk,gridA,gk)
 
 %% (5) 
 B=5000;
 K=zeros(B,20);
for b=1:B 
    Kt=zeros(1:21);
    Kt(1)=int64( find(grid_k>=3,1, 'first') );
    At=zeros(1:21);
    At(1)=int64( find(grid_A>=1,1, 'first') );
    for t=1:20
        At(t+1) = Gen_yind(prob,A(t));
        Kt(t+1) = gk_ind(K(t),A(t));
    end
    K(b,:)=grid_k(Kt(2:20));
end
meanK=mean(K,1);
x=1:20;
plot(x, MeanK, 'b')

%% (6)
% Yes, a stable distribution exists because the technology is not growing. 
% Suppose the population size is N and after T periods the economy reaches
% a stable distribution. 
N=2000;
StableDistr=zeros(N,1); 
T=2000;
for n=1:N
        Kt=zeros(1:T+1);
        Kt(1)=int64( find(grid_k>=3,1, 'first') );
        At=zeros(1:T+1);
        At(1)=int64( find(grid_A>=1,1, 'first') );
        for t=1:T
          At(t+1) = Gen_yind(prob,A(t));
          Kt(t+1) = gk_ind(K(t),A(t));
        end
        StableDistr(n)=Kt(T+1);
end
histogram(StableDistr,40)
    
        

            
             
                 
    













