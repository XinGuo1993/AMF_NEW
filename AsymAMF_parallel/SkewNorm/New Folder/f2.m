function f2 = f2(x,A,B,C)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% function for J_log
    f2 =  exp(B*x - C*x.^2 + A*log(x) + log(log(x)));
end


