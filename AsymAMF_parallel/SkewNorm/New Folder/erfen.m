A = 50;
B = 0;
C = 50;
f=@(x)1./(x.*log(x))+A./x+B-2*C.*x;
a = 0.001;
b = 1-0.001;
c = 1+0.001;
d = 10;
threshold = 1e-6;
%确定 初始点
fa=f(a);
fb=f(b);
fc=f(c);
fd=f(d);
while fa<0
    a = a/2
end
while fb>0
    b = (1+b)/2
end
while fc<0
    c = (1+c)/2
end
while fd>0
    d = 2*d
end
while min(abs(fa),abs(fb))>threshold
    tem = (a+b)*0.5;
    f_tem = f(tem);
    if f_tem>0
        a = tem
        fa = f_tem;
    else
        b = tem
        fb = f_tem; 
    end
    
end

while min(abs(fc),abs(fd))>threshold
    tem = (c+d)*0.5;
    f_tem = f(tem);
    if f_tem>0
        c = tem
        fc = f_tem;
    else
        d = tem
        fd = f_tem;
    end
end

if fa < -fb
    tem1 = a;
else
    tem1 = b;
end
if fc < -fd
    tem2 = c;
else
    tem2 = d;
end
