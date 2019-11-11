%% Question 3 
%% (1) 
% The problem is 
%
% $$ \max_{C_t, K_t, N_t} \beta^t log C_t $$ subject to 
% $$ C_t +K_{t+1} = A_t ^{1-\alpha} K_t^\alpha N_t^{1-\alpha} $$; 
% $$ A_t = (1+g)^t A_0 $$; 
% $$ C_t \geq 0 $$; $$ K_{t+1} \geq 0 $$; $$ N_t \in [0, \bar{N} ] $$. 
% 
% Euler equation is 
%
% $$ \frac{1}{C_t} = \frac{\beta \big( \alpha A_t ^{1-\alpha} K_t^{\alpha-1} \bar{N} ^{1-\alpha} +1 -\delta ) } { C_{t+1} } $$.
% 
% Resource Constraint is 
%
% $$ C_t +K_{t+1} = A_t ^{1-\alpha} K_t^\alpha \bar{N}^{1-\alpha} $$. 
% 
% Define that $$ k_t = K_t (1+g)^{-t} $$ and $$ c_t = C_t (1+g)^{-t} $$ for
% all t. Then 
%
% $$ \frac{ (1+g) c_{t+1} } { c_t} = \beta \big( \alpha A_0^{1-\alpha}
% k_{t+1}^{\alpha -1} \bar{N}^{1-\alpha } + 1 - \delta \big) $$;
%
% $$ (1+g) k_{t+1} =  A_0^{1-\alpha} k_{t}^{\alpha} \bar{N}^{1-\alpha } +
% (1-\delta) k_t -c_t $$. 
%
% Assuming $$ 1+g > \beta (1-\delta ) $$, these equations give a unique
% balanced path, in which $$ k_{t+1} = k_t = k^* >0 $$ and $$ c_{t+1} = c_t
% = c^* >0 $$, with 
%
% $$ k^* = \Big( \frac{ \alpha \beta } { 1+g -\beta (1-\delta)}
% \Big)^\frac{1}{1-\alpha} A_0 \bar{N} $$;
%
% $$ c^* = \Bigg( \big( \frac{ \alpha \beta } { 1+g -\beta (1-\delta)}
% \big)^\frac{\alpha}{1-\alpha} - (\delta +g) \big( \frac{ \alpha \beta } { 1+g -\beta (1-\delta)}
% \big)^\frac{1}{1-\alpha} \Bigg) A_0 \bar{N} $$.
%
%% (2) 
% The phase diagram is like following: 
% 
% <<phase.PNG>>
% 
% The saddle path has the form: 
%
% <<saddle.PNG>>
%
clear
g=0.2;
alpha= 0.3;
beta =0.9; 
delta = 0.5;
k0=1;
barN=1;
A0= 1;
ke = ( alpha*beta/ ( 1+g -beta*(1-delta)))^(1/(1-alpha)) * A0 * barN;
ce = ( ( alpha*beta/ ( 1+g -beta*(1-delta)))^(alpha/(1-alpha)) - ... 
(delta +g)*( alpha*beta/ ( 1+g -beta*(1-delta)))^(1/(1-alpha)) )* A0 * barN;
c0=0.48:0.0001:5;
B=size(c0,2);
Accept=zeros(B,1);
T=6;
k=zeros(T,1);
c=zeros(T,1);
k(1)=k0;
for b=1:B
    c(1)=c0(b);
    for t=1:T
        k(t+1)= (A0^(1-alpha)*k(t)^alpha*barN^(1-alpha) +(1-delta)*k(t)-c(t))/(1+g);
        c(t+1)= beta*( alpha* A0^(1-alpha)*k(t)^(alpha-1)*barN^(1-alpha)+1-delta)*c(t)/(1+g);
        if k(t+1)<ke && c(t+1)>ce 
            break;
        else if k(t+1)>ke && c(t+1) <ce
                break;
            end
        end
    end
    if t==T
        Accept(b)=1;
    end
end
lower= find(Accept,1,'first');
upper= find(Accept,1,'last');
C0 = (c0(lower)+c0(upper))/2;
fprintf('The approximate consumption on saddle path is %1.4f', C0);

%% 
% Comment: You need to set T to be small. A large T will eliminate every
% possible C0. 







