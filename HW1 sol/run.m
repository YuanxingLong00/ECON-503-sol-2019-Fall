%% (2) 
clear
gamma=1;
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
err= 1e-4;
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
    i
end
save('EGM','grid_A','grid_k','gc','gk','gk_ind','prob','invdist','Na','Nk');
