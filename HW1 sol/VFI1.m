R=1.06;
rho=0.8;
Ny=15;
[prob,grid_y,invdist]=tauchen(Ny,1.2,rho,1); % Approximate AR(1) by Markov Process
grid_a= 0.25:0.25:100;
Na= size(grid_a,2); 
V_ini=zeros(Na,Ny);% initial guess for V
V=zeros(Na,Ny);
pol=zeros(Na,Ny);
polI=zeros(Na,Ny);
err= 10^(-6);
iter=1000;
for i=1:iter
     for k=1:Ny
          temp_y= beta* V_ini*(prob(k,:)');
          for j=1:Na
              temp1= log( grid_y(k) + R* grid_a(j) - grid_a');
              temp1( imag(temp1) ~=0) = -Inf;
              temp = temp_y +temp1;
              [M,I]=max(temp,[],1);
              polI(j,k)=I;
              pol(j,k) =grid_a(I);
              V(j,k) =M;
          end
     end
     if norm(V-V_ini)< err
         fprintf('Number of iteratios: %d, Convergence Achieved. \n', i);
         break
     else
         V_ini=V; % update the g=previous guess
     end
end
                       