function ynew_i = Gen_yind(prob, yold_i)
yold_i = int64(  yold_i );
trans= prob(yold_i,:);
u=rand;
CDF=cumsum(trans,2);
ynew_i= int64( find(CDF>=u,1,'first') );



