%% main program
clear
load saxs_data.mat
load saxs_variables.mat
n = numel(saxs_data)/2;
Rg_mod = [];
n_glim_ls = [];
for i = 1:1:n
    s = saxs_data{1,i}(:,1);
    i_s = saxs_data{1,i}(:,2);
    %% find number of points below Guinier limit
    islim = 0;
    n_glim = 1;
    glim = 1.3/Rg(i);
    while islim == 0
        if s(n_glim) < glim
            n_glim = n_glim + 1;
        else
            islim = 1;
        end
    end
    %% find intercept of fitting line in Guinier plot
    n_linear = int16(0.2*n_glim); % exclude high noise region
    guifit = fit(s(n_linear:n_glim).^2, log(i_s(n_linear:n_glim)), 'poly1');
    i0 = guifit.p2;
    i0 = exp(i0);
    Rg_mod(i) = sqrt(guifit.p1*-3);
    %% update n_glim according to newly refined Rg
    islim = 0;
    n_glim = 1;
    glim = 1.3/Rg_mod(i);
    while islim == 0
        if s(n_glim) < glim
            n_glim = n_glim + 1;
        else
            islim = 1;
        end
    end
    %% normalize data
    i_s = i_s./i0;
    %% model high angle part
    s_ha = s(n_glim:end);
    i_s_ha = i_s(n_glim:end);
    ft = fittype('A*x^-4 + B*x^-3 + C*x^-2 + D*x^-1 + E');
    hafit = fit(s_ha, i_s_ha, ft);
    A(i) = hafit.A;
    B(i) = hafit.B;
    C(i) = hafit.C;
    D(i) = hafit.D;
    E(i) = hafit.E;
    s_jnt = s(n_glim);
    i_jnt_g = exp(-1/3*Rg_mod(i)^2*s_jnt^2);
    i_jnt_ha = A(i)*s_jnt^-4 + B(i)*s_jnt^-3 + C(i)*s_jnt^-2 + D(i)*s_jnt^-1 + E(i);
    contrat = i_jnt_g/i_jnt_ha;
    A(i) = contrat*A(i);
    B(i) = contrat*B(i);
    C(i) = contrat*C(i);
    D(i) = contrat*D(i);
    E(i) = contrat*E(i);
    
    %% generate data-fit plot
%     plot(s, i_s)
%     hold on
%     plot(s(1:n_glim), exp(-1/3*Rg_mod(i)^2.*s(1:n_glim).^2))
%     plot(s_ha, A(i)*s_ha.^-4 + B(i)*s_ha.^-3 + C(i)*s_ha.^-2 + D(i)*s_ha.^-1 + E(i))
%     plot(s_jnt, i_jnt_g, 'o');
%     hold off
%     axis([0, 5, 0, 1])
%     title([num2str(i), ': ', saxs_data{2,i}]);
%     pause();
end
%% save all fitting curves
saxs_res = zeros(n, 5001);
s = 0:0.001:5;
for i = 1:1:n
    islim = 0;
    n_glim = 1;
    glim = 1.3/Rg_mod(i);
    while islim == 0
        if s(n_glim) < glim
            n_glim = n_glim + 1;
        else
            islim = 1;
        end
    end
    saxs_res(i, 1:n_glim) = exp(-(1/3)*Rg_mod(i)^2*s(1:n_glim).^2);
    s_ha = s(n_glim+1:end);
    saxs_res(i, n_glim+1:end) = A(i)*s_ha.^-4 + B(i)*s_ha.^-3 + C(i)*s_ha.^-2 + D(i)*s_ha.^-1 + E(i);
end
%% plot all fits in one graph
for i = 1:n
    loglog(s, saxs_res(i,:));
    hold on
end
%% average all fitting curves
fit_ave = zeros(1, 5001);
for i = 1:1:5001
    fit_ave(i) = mean(saxs_res(1:end, i));
end
plot(s, fit_ave);
%% fit averaged curve
glim = 1.3/mean(Rg_mod);
islim = 0;
n_glim = 1;
while islim == 0
    if s(n_glim) < glim
        n_glim = n_glim + 1;
    else
        islim = 1;
    end
end
n_linear = int16(0.2*n_glim); 
s = s'; % transpose for fit
fit_ave = fit_ave'; % transpose for fit
guifit = fit(s(n_linear:n_glim).^2, log(i_s(n_linear:n_glim)), 'poly1');
i0 = guifit.p2;
i0 = exp(i0);
Rg_final = sqrt(guifit.p1*-3);
%% update n_glim
glim = 1.3/Rg_final;
islim = 0;
n_glim = 1;
while islim == 0
    if s(n_glim) < glim
        n_glim = n_glim + 1;
    else
        islim = 1;
    end
end
%% normalize data
i_s = i_s./i0;
%% model high angle part 
s_ha = s(n_glim:end);
i_s_ha = fit_ave(n_glim:end);
ft = fittype('A*x^-4 + B*x^-3 + C*x^-2 + D*x^-1 + E');
hafit = fit(s_ha, i_s_ha, ft);
A_final = hafit.A;
B_final = hafit.B;
C_final = hafit.C;
D_final = hafit.D;
E_final = hafit.E;
s_jnt = s(n_glim);
i_jnt_g = exp(-1/3*Rg_final^2*s_jnt^2);
i_jnt_ha = A_final*s_jnt^-4 + B_final*s_jnt^-3 + C_final*s_jnt^-2 + D_final*s_jnt^-1 + E_final;
contrat = i_jnt_g/i_jnt_ha;
A_final = contrat*A_final;
B_final = contrat*B_final;
C_final = contrat*C_final;
D_final = contrat*D_final;
E_final = contrat*E_final;
%% generate final fit plot
figure();
loglog(s(1:n_glim), exp(-1/3*Rg_final^2.*s(1:n_glim).^2))
hold on
loglog(s_ha, A_final*s_ha.^-4 + B_final*s_ha.^-3 + C_final*s_ha.^-2 + D_final*s_ha.^-1 + E_final)
hold off
axis([0, 5, 0, 1])
%% output parameters
disp([Rg_final, A_final, B_final, C_final, D_final, E_final]);

%% MW histogram
figure();
histogram(mw,'BinWidth',50);
