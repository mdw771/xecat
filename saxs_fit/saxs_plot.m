AA = mean(A);
BB = mean(B);
CC = mean(C);
DD = mean(D);
EE = mean(E);
RR = mean(Rg);
x = 0:0.01:5;
islim = 0;
n_glim = 1;
glim = 1.3/RR;
while islim == 0
    if x(n_glim) < glim
        n_glim = n_glim + 1;
    else
        islim = 1;
    end
end
x_ha = x(n_glim:end);
plot(x, exp(-(1/3)*RR^2*x.^2));
hold on
plot(x_ha, AA*x_ha.^-4 + BB*x_ha.^-3 + CC*x_ha.^-2 + DD*x_ha.^-1 + EE);
line([1.3/RR, 1.3/RR], [0, 1]);
hold off