function sol= runnested(knew, a, prob, gc, grid_A,grid_k,X0)

x=X0;
tol=1e-6;
maxit=500;
for i=1:maxit
    fx=Obj(x);
    x = x - Jac(x)\fx;
    if norm(fx)< tol
        break
    end
end
sol=x;



function F=Obj(unknowns)
c=unknowns(1);
k=unknowns(2);
gamma=0.5;
alpha=0.3;
beta=0.9;
delta=0.95;
F=zeros(2,1);
F(1)= 1-c^gamma*  beta* ( ( alpha* grid_A.^(alpha-1).* grid_k(knew)^(alpha-1)...
    + 1-delta ) ./( gc(knew,:).^gamma) ) *prob(a,:)';
F(2) = grid_k(knew)+c-grid_A(a)*k^alpha - (1-delta)*k;
end

    
    

function J = Jac(unknowns)
c=unknowns(1);
k=unknowns(2);
gamma=0.5;
alpha=0.3;
beta=0.9;
delta=0.95;
J = ones(2,2);
J(1,1)= gamma* c^(gamma-1) * beta* ( ( alpha* grid_A.^(alpha-1).* grid_k(knew)^(alpha-1)...
    + 1-delta ) ./( gc(knew,:).^gamma) ) *prob(a,:)';
J(1,2) =0;
J(2,1) = 1;
J(2,2) =-grid_A(a)*alpha* k^(alpha-1)- (1-delta);
end

end





