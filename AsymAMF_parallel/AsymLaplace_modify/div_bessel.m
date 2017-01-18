function result = div_bessel(mu,Z,h)
%近似求解bessel比例
if nargin~=3
    h=1;
end
mu0 = mu;
flag = 0;
if (besselk(mu+h,Z) == Inf) | (besselk(mu,Z) == Inf)
  flag = 1;
end

while flag ==1
    mu = mu*0.5;
    if (besselk(mu+h,Z) ~= Inf) && (besselk(mu,Z) ~= Inf)
        flag = 2;
    end
end

if flag ==0
    result = besselk(mu+h,Z)./besselk(mu,Z);
elseif flag == 2
    x = linspace(mu*0.5,mu);
    Z = ones(1,100)*Z;
    y = besselk(mu+h,Z)./besselk(mu,Z);
    p = polyfit(x,y,1);
    result=p(1)*mu0+p(2);    
end
end