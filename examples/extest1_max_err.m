function e = extest1_max_err(p)
A = 1.0;
alpha = 2.0;
C0 = exp(4*alpha) - 1;

e = 0;
for i = 1:length(p.T)
    t = p.T(i);
    u = 4*alpha*exp(alpha*t) / (A*C0);
    e = max(e, abs(u-p.U(1,i)));
end

end
