clearvars;%close all;

% Preparing the data
dataGRdep = readtable('EuropeanGradientData.xlsx','Sheet','EuropeanGradientData');
data2 = readtable('EuropeanGradientData.xlsx','Sheet','EG_SitesDataset');

micGrow = dataGRdep.b_c + dataGRdep.f_c; % bacteria and fungi combined
micResp = dataGRdep.resp; 
tempK = dataGRdep.temp + 273.15;

% Preparing variables
nSites = 70;
nModels = 4;
nModels2 = 5;

W_gr = zeros(nSites,7);   
H_gr = zeros(nSites,5);   
M_gr = zeros(nSites,4);   
D_gr = zeros(nSites,6);

W_re = zeros(nSites,8);
H_re = zeros(nSites,6);
M_re = zeros(nSites,5);
D_re = zeros(nSites,8);

% One-step calibration (MMRTD_comb)
D_co = zeros(nSites,8);

R2_gr   = zeros(nSites,nModels);
RMSE_gr = zeros(nSites,nModels);
AIC_gr  = zeros(nSites,nModels);

R2_re   = zeros(nSites,nModels);
RMSE_re = zeros(nSites,nModels);
AIC_re  = zeros(nSites,nModels);

R2_co   = zeros(nSites,1);
RMSE_co = zeros(nSites,1);
AIC_co  = zeros(nSites,1);

R2_CUE = zeros(nSites,nModels2);
RMSE_CUE = zeros(nSites,nModels2);
AIC_CUE = zeros(nSites,nModels2);

MAST = zeros(70,1);

% Preparing for loop
posTdep = 1;
posSite = 1;
T_fit = linspace(273.15 - 5, 273.15 + 55, 100);
R = 8.314462618; Kb = 1.380649e-23;h = 6.62607015e-34;

for ii = 1
    if ii == 3
        posTdep = posTdep+10;
    end
    if ii == 25
        posTdep = posTdep+10;
        posSite = posSite+1;
    end

    SoilID = posTdep;

    MAST(ii) = data2(posSite,14).Variables;

    % Preparing data
    T_ = tempK(posTdep : posTdep+9);
    G_ = micGrow(posTdep : posTdep+9);
    R_ = micResp(posTdep: posTdep+9);
    
    nanIn_gr = isnan(G_); % identify NaN
    growth = G_(~nanIn_gr); % filtered growth data
    temp_gr = T_(~nanIn_gr);
    
    nanIn_re = isnan(R_);
    resp = R_(~nanIn_re); 
    temp_re = T_(~nanIn_re); 
    
    gr_re = [growth; resp];
    temp_regr = [temp_re; temp_gr];

    nanIn_comb = isnan(G_)|isnan(R_);
    Gcomb= G_(~nanIn_comb);
    Rcomb= R_(~nanIn_comb);
    temp_comb = T_(~nanIn_comb);

    % Initial values - Growth
    ini_Huang = [9.3 2.6e3 13 7e-5 318];
    ini_Walker = [5.7e4 -70 120 -5.4e5 1.8e5 360 180];
    ini_MMRT = [3.7e6 -3.2e3 863 2.85e03];

    % Function - Growth 
    funHuang = @(T,P) P(1)*T.*exp(-(P(2)/R./T).^P(3)).*(1-exp(P(4)*(T-P(5))));
    funCpT = @(T,P) (P(3)+P(4)*exp(-P(5)*(1-T/P(6))/R./T)) ./ (1+exp(-P(5)*(1-T/P(6))/R./T));
    funWalker = @(T,P) (Kb/h*T) .* exp( -(P(1) - T*P(2) + funCpT(T,P).*(T-P(7)-T.*log(T/P(7))) ) /R./T);
    funMMRT = @(T,P) Kb/h.*T.*exp((-P(1)+P(2).*(T-P(3)))./(R.*T)+(P(4)+P(2).*(log(T)-log(P(3))))/R);
    funGrow_MMRTD = @(T,P) (Kb*T/h) .* (exp((-P(1)+P(2).*(T-P(3)))./(R.*T)+(P(4)+P(2).*(log(T)-log(P(3))))/R) - exp((-P(5)./(R.*T))+(P(6)/R)));

    % Parameter estimation - Growth
    param_MMRT_aux = Regression_Tmodel_s1(temp_gr,growth,funMMRT,ini_MMRT);
    param_MMRT = Regression_Tmodel_s2(temp_gr,growth,funMMRT,param_MMRT_aux);
    param_Huang = Regression_Tmodel_s2(temp_gr,growth,funHuang,ini_Huang);
    param_Walker = Regression_Tmodel_s2(temp_gr, growth, funWalker, ini_Walker);
    param_MMRTD = Regression_Tmodel_s2(temp_gr, growth, funGrow_MMRTD, [param_MMRT_aux  6.9e4 -32]); 

    % Plotting - Growth
    figure(ii)
    subplot(1,3,1);plot(temp_gr, growth, 'ok','Displayname','Observed Data','MarkerSize',10); 
                   hold on
                   plot(T_fit,funMMRT(T_fit,param_MMRT),'g-','DisplayName','MMRT');
                   hold on
                   plot(T_fit,funHuang(T_fit,param_Huang),'y-', 'DisplayName', 'Huang');
                   hold on
                   plot(T_fit,funWalker(T_fit,param_Walker),'c-', 'DisplayName', 'Walker');
                   hold on
                   plot(T_fit,funGrow_MMRTD(T_fit,param_MMRTD),'m-', 'DisplayName','MMRT+D_g');
                   lgd = legend('show');
                   lgd.FontSize = 10; 
                   ylim([0 max(growth)*1.3]);
                   xlabel('Temperature (K)','FontSize', 12)
                   ylabel('Microbial growth','FontSize',12)

    % Initial values - Respiration
    pCF1 = resp(6)/growth(6);
    pCF2 = max(resp)/max(growth);
    iniResp_Walker = [param_Walker pCF1];
    iniResp_Huang = [param_Huang pCF1];
    iniResp_MMRT = [param_MMRT_aux pCF1];
    iniResp_MMRTD = [param_MMRTD(1:6) pCF1 pCF2];
    
    % Function - Respiration
    funResp_Walker = @(T,P) funWalker(T,P)*P(8);
    funResp_Huang = @(T,P) funHuang(T,P)*P(6);
    funResp_MMRT = @(T,P) funMMRT(T,P)*P(5);
    funResp_MMRTD = @(T,P) Kb*T/h.*(exp((-P(1)+P(2).*(T-P(3)))./(R.*T)+(P(4)+P(2).*(log(T)-log(P(3))))/R)*P(7)+ exp((-P(5)./(R.*T))+(P(6)/R))*P(8));

    % Parameter estimation - Respiration
    paramRe_MMRT = Regression_Tmodel_s2(temp_re,resp,funResp_MMRT,iniResp_MMRT);
    paramRe_Huang = Regression_Tmodel_s2(temp_re,resp,funResp_Huang,iniResp_Huang);
    paramRe_Walker = Regression_Tmodel_s2(temp_re, resp, funResp_Walker, iniResp_Walker);
    paramRe_MMRTD = Regression_Tmodel_s2(temp_re, resp, funResp_MMRTD,iniResp_MMRTD);      
    
    % Parameter estimation - MMRTD_comb
    paramComb_MMRTD = Regression_Tmodel_comb(temp_gr,growth,funGrow_MMRTD,temp_re,resp,funResp_MMRTD,iniResp_MMRTD); 

    % Plotting - Respiration
    subplot(1,3,2);plot(temp_re, resp, 'ok','Displayname','Observed Data', 'MarkerSize',10); % observed microbial respiration for the unique soil
                   hold on
                   plot(T_fit,funResp_MMRT(T_fit,paramRe_MMRT),'g-','DisplayName','MMRT');
                   hold on
                   plot(T_fit,funResp_Huang(T_fit,paramRe_Huang),'y-', 'DisplayName', 'Huang');
                   hold on
                   plot(T_fit,funResp_Walker(T_fit,paramRe_Walker),'c-', 'DisplayName', 'Walker');
                   hold on
                   plot(T_fit,funResp_MMRTD(T_fit,paramRe_MMRTD),'m-', 'DisplayName','MMRT+D_r');
                   hold on
                   plot(T_fit,funResp_MMRTD(T_fit,paramComb_MMRTD),'k:','DisplayName','MMRT+Dcomb'); % MMRTD one-step
                   lgd = legend('show');
                   lgd.FontSize = 10;
                   ylim([0 max(resp)*1.5]);
                   xlabel('Temperature (K)','FontSize',12)
                   ylabel('Microbial respiration','FontSize',12)
    subplot(1,3,1); plot(T_fit,funGrow_MMRTD(T_fit,paramComb_MMRTD(1:6)),'k:','DisplayName','MMRT+Dcomb'); % MMRTD one-step      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting CUE
    subplot(1,3,3); plot(temp_comb,(Gcomb./(Gcomb+Rcomb)),'ok','DisplayName','CUE Data','MarkerSize',10);
                    hold on
                    plot(T_fit,funWalker(T_fit,param_Walker)./(funResp_Walker(T_fit,paramRe_Walker)+funWalker(T_fit,param_Walker)),'c-','DisplayName','CUE Walker');
                    hold on
                    plot(T_fit,funHuang(T_fit,param_Huang)./(funResp_Huang(T_fit,paramRe_Huang)+funHuang(T_fit,param_Huang)),'y-','DisplayName','CUE Huang');
                    hold on
                    plot(T_fit,funMMRT(T_fit,param_MMRT)./(funResp_MMRT(T_fit,paramRe_MMRT)+funMMRT(T_fit,param_MMRT)),'g-','DisplayName','CUE MMRT');
                    hold on
                    plot(T_fit,funGrow_MMRTD(T_fit,param_MMRTD)./(funResp_MMRTD(T_fit,paramComb_MMRTD)+funGrow_MMRTD(T_fit,param_MMRTD)),'m-','DisplayName','CUE MMRT_D');
                    hold on
                    plot(T_fit,funGrow_MMRTD(T_fit,paramComb_MMRTD(1:6))./(funResp_MMRTD(T_fit,paramComb_MMRTD)+funGrow_MMRTD(T_fit,paramComb_MMRTD(1:6))),'k:','DisplayName','CUE MMRT_D Combined');
                    xlabel('Temperature (K)','Fontsize',12);
                    ylabel('CUE','FontSize',12);
                    lgd = legend('show');
                    lgd.FontSize = 10;

    % Parameter storing
    W_gr(ii,:) = param_Walker;   
    H_gr(ii,:) = param_Huang;    
    M_gr(ii,:) = param_MMRT;     
    D_gr(ii,:) = param_MMRTD;    % MMRTD two-step   

    W_re(ii,:) = paramRe_Walker;
    H_re(ii,:) = paramRe_Huang;
    M_re(ii,:) = paramRe_MMRT;
    D_re(ii,:) = paramRe_MMRTD; % MMRTD two-step     
    
    D_co(ii,:) = paramComb_MMRTD;  % MMRTD_comb one-step
    
    maxG = max(growth);
    maxR = max(resp);
    [R2_gr(ii,1), RMSE_gr(ii,1), AIC_gr(ii,1)] = funModEval(growth/maxG, funWalker(temp_gr,param_Walker)/maxG,7);            
    [R2_gr(ii,2), RMSE_gr(ii,2), AIC_gr(ii,2)] = funModEval(growth/maxG, funHuang(temp_gr,param_Huang)/maxG,5);            
    [R2_gr(ii,3), RMSE_gr(ii,3), AIC_gr(ii,3)] = funModEval(growth/maxG, funMMRT(temp_gr,param_MMRT)/maxG,4);            
    [R2_gr(ii,4), RMSE_gr(ii,4), AIC_gr(ii,4)] = funModEval(growth/maxG, funGrow_MMRTD(temp_gr,param_MMRTD)/maxG,6);            
    [R2_gr(ii,5), RMSE_gr(ii,5), AIC_gr(ii,5)] = funModEval(growth/maxG, funGrow_MMRTD(temp_gr,paramComb_MMRTD)/maxG,6);            
    
    [R2_re(ii,1), RMSE_re(ii,1), AIC_re(ii,1)] = funModEval(resp/maxR, funResp_Walker(temp_re,paramRe_Walker)/maxR,8);            
    [R2_re(ii,2), RMSE_re(ii,2), AIC_re(ii,2)] = funModEval(resp/maxR, funResp_Huang(temp_re,paramRe_Huang)/maxR,6);            
    [R2_re(ii,3), RMSE_re(ii,3), AIC_re(ii,3)] = funModEval(resp/maxR, funResp_MMRT(temp_re,paramRe_MMRT)/maxR,5);            
    [R2_re(ii,4), RMSE_re(ii,4), AIC_re(ii,4)] = funModEval(resp/maxR, funResp_MMRTD(temp_re,paramRe_MMRTD)/maxR,8);            
    [R2_re(ii,5), RMSE_re(ii,5), AIC_re(ii,5)] = funModEval(resp/maxR, funResp_MMRTD(temp_re,paramComb_MMRTD)/maxR,8);            
    
    [R2_co(ii,1), RMSE_co(ii,1), AIC_co(ii,1)] = funModEval([growth/maxG;resp/maxR], [funWalker(temp_gr,param_Walker)/maxG;funResp_Walker(temp_re,paramRe_Walker)/maxR],15);      
    [R2_co(ii,2), RMSE_co(ii,2), AIC_co(ii,2)] = funModEval([growth/maxG;resp/maxR], [funHuang(temp_gr,param_Huang)/maxG;funResp_Huang(temp_re,paramRe_Huang)/maxR],11);      
    [R2_co(ii,3), RMSE_co(ii,3), AIC_co(ii,3)] = funModEval([growth/maxG;resp/maxR], [funMMRT(temp_gr,param_MMRT)/maxG;funResp_MMRT(temp_re,paramRe_MMRT)/maxR],9);      
    [R2_co(ii,4), RMSE_co(ii,4), AIC_co(ii,4)] = funModEval([growth/maxG;resp/maxR], [funGrow_MMRTD(temp_gr,param_MMRTD)/maxG;funResp_MMRTD(temp_re,paramRe_MMRTD)/maxR],14);      
    [R2_co(ii,5), RMSE_co(ii,5), AIC_co(ii,5)] = funModEval([growth/maxG;resp/maxR], [funGrow_MMRTD(temp_gr,paramComb_MMRTD)/maxG;funResp_MMRTD(temp_re,paramComb_MMRTD)/maxR],8);      

    [R2_CUE(ii,1), RMSE_CUE(ii,1), AIC_CUE(ii,1)] = funModEval((Gcomb./(Gcomb+Rcomb)), funWalker(temp_comb,param_Walker)./(funResp_Walker(temp_comb,paramRe_Walker)+funWalker(temp_comb,param_Walker)),15);      
    [R2_CUE(ii,2), RMSE_CUE(ii,2), AIC_CUE(ii,2)] = funModEval((Gcomb./(Gcomb+Rcomb)), funHuang(temp_comb,param_Huang)./(funResp_Huang(temp_comb,paramRe_Huang)+funHuang(temp_comb,param_Huang)),11);      
    [R2_CUE(ii,3), RMSE_CUE(ii,3), AIC_CUE(ii,3)] = funModEval((Gcomb./(Gcomb+Rcomb)), funMMRT(temp_comb,param_MMRT)./(funResp_MMRT(temp_comb,paramRe_MMRT)+funMMRT(temp_comb,param_MMRT)),9);      
    [R2_CUE(ii,4), RMSE_CUE(ii,4), AIC_CUE(ii,4)] = funModEval((Gcomb./(Gcomb+Rcomb)), funGrow_MMRTD(temp_comb,param_MMRTD)./(funResp_MMRTD(temp_comb,paramComb_MMRTD)+funGrow_MMRTD(temp_comb,param_MMRTD)),14);      
    [R2_CUE(ii,5), RMSE_CUE(ii,5), AIC_CUE(ii,5)] = funModEval((Gcomb./(Gcomb+Rcomb)), funGrow_MMRTD(temp_comb,paramComb_MMRTD(1:6))./(funResp_MMRTD(temp_comb,paramComb_MMRTD)+funGrow_MMRTD(temp_comb,paramComb_MMRTD(1:6))),8);     

    posTdep = posTdep+10; % New iteration
    posSite = posSite + 1;
end

% Kruskal-Wallis & Dunn-Sidak
modelNames = {'Walker','Huang','MMRT','MMRT_D','MMRT_Dcomb'};
metrics = {'AIC','R2','RMSE'};
varNames = {'Model_i','Model_j','Lower_CI','Mean_Diff','Upper_CI','p_value'};

% Growth
data_growth = {AIC_gr, R2_gr, RMSE_gr};

for i = 1:numel(metrics)
    [~, ~, stats] = kruskalwallis(data_growth{i}, modelNames, 'off');
    c = multcompare(stats, 'CType', 'dunn-sidak');
    T = array2table(c, 'VariableNames', varNames);
    T.Model_i = string(modelNames(T.Model_i))';
    T.Model_j = string(modelNames(T.Model_j))';
    results_growth.(metrics{i}) = T;
end

%Respiration
data_resp = {AIC_re, R2_re, RMSE_re};

for i = 1:numel(metrics)
    [~, ~, stats] = kruskalwallis(data_resp{i}, modelNames, 'off');
    c = multcompare(stats, 'CType', 'dunn-sidak');
    T = array2table(c, 'VariableNames', varNames);
    T.Model_i = string(modelNames(T.Model_i))';
    T.Model_j = string(modelNames(T.Model_j))';
    results_resp.(metrics{i}) = T;
end

% MMRT_Dcomb
data_co = {AIC_co, R2_co, RMSE_co};

for i = 1:numel(metrics)
    [~, ~, stats] = kruskalwallis(data_co{i}, modelNames, 'off');
    c = multcompare(stats, 'CType', 'dunn-sidak');
    T = array2table(c, 'VariableNames', varNames);
    T.Model_i = string(modelNames(T.Model_i))';
    T.Model_j = string(modelNames(T.Model_j))';
    results_co.(metrics{i}) = T;
end
    
% CUE
data_CUE = {AIC_CUE, R2_CUE, RMSE_CUE};

for i = 1:numel(metrics)
    [~, ~, stats] = kruskalwallis(data_CUE{i}, modelNames, 'off');
    c = multcompare(stats, 'CType', 'dunn-sidak');
    T = array2table(c, 'VariableNames', varNames);
    T.Model_i = string(modelNames(T.Model_i))';
    T.Model_j = string(modelNames(T.Model_j))';
    results_CUE.(metrics{i}) = T;
end

% Spearman Order: [Rho P value]
Walker_Spearman_gr  = spearmanMatrix(MAST, W_gr);
Huang_Spearman_gr   = spearmanMatrix(MAST, H_gr);
MMRT_Spearman_gr    = spearmanMatrix(MAST, M_gr);
MMRTD_Spearman_gr   = spearmanMatrix(MAST, D_gr);

Walker_Spearman_re  = spearmanMatrix(MAST, W_re);
Huang_Spearman_re   = spearmanMatrix(MAST, H_re);
MMRT_Spearman_re    = spearmanMatrix(MAST, M_re);
MMRTD_Spearman_re   = spearmanMatrix(MAST, D_re);

MMRTD_Spearman_co  = spearmanMatrix(MAST, D_co);

% Pattern score
Walker_Pat_gr = -log(Walker_Spearman_gr(:,2)) .* (Walker_Spearman_gr(:,1)).^2;
Huang_Pat_gr  = -log(Huang_Spearman_gr(:,2)) .* (Huang_Spearman_gr(:,1)).^2;
MMRT_Pat_gr   = -log(MMRT_Spearman_gr(:,2)) .* (MMRT_Spearman_gr(:,1)).^2;
MMRTD_Pat_gr  = -log(MMRTD_Spearman_gr(:,2)) .* (MMRTD_Spearman_gr(:,1)).^2;

Walker_Pat_re = -log(Walker_Spearman_re(:,2)) .* (Walker_Spearman_re(:,1)).^2;
Huang_Pat_re  = -log(Huang_Spearman_re(:,2)) .* (Huang_Spearman_re(:,1)).^2;
MMRT_Pat_re   = -log(MMRT_Spearman_re(:,2)) .* (MMRT_Spearman_re(:,1)).^2;
MMRTD_Pat_re  = -log(MMRTD_Spearman_re(:,2)) .* (MMRTD_Spearman_re(:,1)).^2;

MMRTD_Pat_co  = -log(MMRTD_Spearman_co(:,2)) .* (MMRTD_Spearman_co(:,1)).^2;

% Fitting linear regression
% Growth
for c = 1:7
    mdl = fitlm(MAST, W_gr(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);
    pval_slope_W_gr (c) = mdl.Coefficients.pValue(2); % slope p-value

    p_W_gr(:,c) = [slope; intercept];
    r_W_gr(:,c) = slope * MAST + intercept;
end

for c = 1:5
    mdl = fitlm(MAST, H_gr(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_H_gr(:,c) = [slope; intercept];
    r_H_gr(:,c) = slope * MAST + intercept;
end

for c = 1:4
    mdl = fitlm(MAST, M_gr(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_M_gr(:,c) = [slope; intercept];
    r_M_gr(:,c) = slope * MAST + intercept;
end

for c = 1:6
    mdl = fitlm(MAST, D_gr(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_D_gr(:,c) = [slope; intercept];
    r_D_gr(:,c) = slope * MAST + intercept;
end
% Respiration
for c = 1:8
    mdl = fitlm(MAST, W_re(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_W_re(:,c) = [slope; intercept];
    r_W_re(:,c) = slope * MAST + intercept;
end

for c = 1:6
    mdl = fitlm(MAST, H_re(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_H_re(:,c) = [slope; intercept];
    r_H_re(:,c) = slope * MAST + intercept;
end

for c = 1:5
    mdl = fitlm(MAST, M_re(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_M_re(:,c) = [slope; intercept];
    r_M_re(:,c) = slope * MAST + intercept;
end

for c = 1:8
    mdl = fitlm(MAST, D_re(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_D_re(:,c) = [slope; intercept];
    r_D_re(:,c) = slope * MAST + intercept;
end

% MMRT_Dcomb
for c = 1:8
    mdl = fitlm(MAST, D_co(:,c));

    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    p_D_co(:,c) = [slope; intercept];
    r_D_co(:,c) = slope * MAST + intercept;
end

% Pval and Rho storing
Pval_W_gr = Walker_Spearman_gr(:,2);
Rho_W_gr  = Walker_Spearman_gr(:,1);
Pval_H_gr = Huang_Spearman_gr(:,2);
Rho_H_gr  = Huang_Spearman_gr(:,1);
Pval_M_gr = MMRT_Spearman_gr(:,2);
Rho_M_gr  = MMRT_Spearman_gr(:,1);
Pval_D_gr = MMRTD_Spearman_gr(:,2);
Rho_D_gr  = MMRTD_Spearman_gr(:,1);

Pval_W_re = Walker_Spearman_re(:,2);
Rho_W_re  = Walker_Spearman_re(:,1);
Pval_H_re = Huang_Spearman_re(:,2);
Rho_H_re  = Huang_Spearman_re(:,1);
Pval_M_re = MMRT_Spearman_re(:,2);
Rho_M_re  = MMRT_Spearman_re(:,1);
Pval_D_re = MMRTD_Spearman_re(:,2);
Rho_D_re  = MMRTD_Spearman_re(:,1);

Pval_D_co = MMRTD_Spearman_co(:,2);
Rho_D_co  = MMRTD_Spearman_co(:,1);

% Plotting
% Growth
% Walker
figure;
    subplot(2,4,1);
    n = 1;
    plot(MAST, r_W_gr(:,n), 'r-')
    hold on
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    ylim([0 1.5e5]);
    title('P(1) vs MAST','FontSize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kJ/mol)','FontSize',10)
    
    subplot(2,4,2);
    n = 2;
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_W_gr(:,n), 'r-');
    ylim([-1000 2500]);
    title('P(2) vs MAST','Fontsize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,4,3);
    n = 3;
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_W_gr(:,n), 'r-');
    ylim([-500 700]);
    title('P(3) vs MAST','Fontsize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)
   
    subplot(2,4,4);
    n = 4;
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_W_gr(:,n), 'r-');
    ylim([-2e6 0.5e6]);
    title('P(4) vs MAST','FontSize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity low (kJ/mol/K)','FontSize',10)

    subplot(2,4,5);
    n = 5;
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_W_gr(:,n), 'r-');
    title('P(5) vs MAST','Fontsize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity high (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,6);
    n = 6;
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_W_gr(:,n), 'r-');
    title('P(6) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy difference (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,7);
    n = 7;
    plot(MAST, W_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_W_gr(:,n), 'r-');
    title('P(7) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Tc (K)','FontSize',10)

    sgtitle('Walker (Growth) Parameters vs. Mean Annual Soil Temperature (MAST)');

% Huang
figure;
    subplot(2,3,1);
    n = 1;
    plot(MAST, H_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_H_gr(:,n), 'y-');
    ylim([-50 300]);
    title('P(1) vs MAST','FontSize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Frequency factor','FontSize',10)
    
    subplot(2,3,2);
    n = 2;
    plot(MAST, H_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_H_gr(:,n), 'y-');
    ylim([2300 3200]);
    title('P(2) vs MAST','FontSize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Energy term (J/mol)','FontSize',10)
    
    subplot(2,3,3);
    n = 3;
    plot(MAST, H_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_H_gr(:,n), 'y-');
    title('P(3) vs MAST','FontSize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting exponent a','FontSize',10)
    
    subplot(2,3,4);
    n = 4;
    plot(MAST, H_gr(:,n), 'ko', 'MarkerSize', 6);
    hold on
    plot(MAST, r_H_gr(:,n), 'y-');
    ylim([-0.05 0.05]);
    title('P(4) vs MAST','FontSize',12);
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Tmax (K)','FontSize',10)
    
    subplot(2,3,5);
    n = 5;
    plot(MAST, H_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_H_gr(:,n), 'y-');
    ylim([314 324]);
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting exponent c','FontSize',10)
 
    sgtitle('Huang (Growth) Parameters vs. Mean Annual Soil Temperature (MAST)');

% MMRT
figure;
    subplot(2,2,1);
    n = 1;
    plot(MAST, M_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_M_gr(:,n), 'c-');
    ylim([-5 1.5e7]);
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,2,2);
    n = 2;
    plot(MAST, M_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_M_gr(:,n), 'c-');
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,2,3);
    n = 3;
    plot(MAST, M_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_M_gr(:,n), 'c-');
    ylim([-100 1500]);
    title('P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('T0 (K)','FontSize',10)
    
    subplot(2,2,4);
    n = 4;
    plot(MAST, M_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_M_gr(:,n), 'c-');
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)
   
    sgtitle('MMRT (Growth) Parameters vs. Mean Annual Soil Temperature (MAST)');

%MMRTD
figure;
    subplot(2,3,1);
    n = 1;
    plot(MAST, D_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_gr(:,n), 'm-');
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,3,2);
    n = 2;
    plot(MAST, D_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_gr(:,n), 'm-');
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,3,3);
    n = 3;
    plot(MAST, D_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_gr(:,n), 'm-');
    ylim([300 1400]);
    title('P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('T0 (K)','FontSize',10)
    
    subplot(2,3,4);
    n = 4;
    plot(MAST, D_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_gr(:,n), 'm-');
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)
    
    subplot(2,3,5);
    n = 5;
    plot(MAST, D_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_gr(:,n), 'm-');
    ylim([5.5e4 8.5e4]);
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Difference in Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,3,6);
    n = 6;
    plot(MAST, D_gr(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_gr(:,n), 'm-');
    ylim([-70 0]);
    title('P(6) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Difference in Entropy (J/mol/K)','FontSize',10)
    
    sgtitle('MMRT_D (Growth) Parameters vs. Mean Annual Soil Temperature (MAST)');

% Respiration
% Walker
figure;
    subplot(2,4,1);
    n = 1;
    plot(MAST, W_re(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([-0.5e6 0.5e6]);
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kJ/mol)','FontSize',10)
    
    subplot(2,4,2);
    n = 2;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([-500 500]);
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
   
    subplot(2,4,3);
    n = 3;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([-1300 1500]);
    title('P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting exponent a','FontSize',10)
    
    subplot(2,4,4);
    n = 4;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([-2e6 2e6]);
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,5);
    n = 5;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([3e4 5e5]);
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity high (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,6);
    n = 6;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([0 3000]);
    title('P(6) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy difference (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,7);
    n = 7;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([10 500]);
    title('P(7) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Tc (K)','FontSize',10)
    
    subplot(2,4,8);
    n = 8;
    plot(MAST, W_re(:,n), 'ko');
    hold on
    plot(MAST, r_W_re(:,n), 'r-');
    ylim([10 500]);
    title('P(8) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter','FontSize',10)
    
    sgtitle('Walker (Respiration) Parameters vs. Mean Annual Soil Temperature (MAST)');

% Huang
figure;
    subplot(2,3,1);
    n = 1;
    plot(MAST, H_re(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_H_re(:,n), 'y-');
    ylim([-1000 5000]);
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Frequency factor','FontSize',10)
    
    subplot(2,3,2);
    n = 2;
    plot(MAST, H_re(:,n), 'ko');
    hold on
    plot(MAST, r_H_re(:,n), 'y-');
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,3,3);
    n = 3;
    plot(MAST, H_re(:,n), 'ko');
    hold on
    plot(MAST, r_H_re(:,n), 'y-');
    title('P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('T0 (K)','FontSize',10)
    
    subplot(2,3,4);
    n = 4;
    plot(MAST, H_re(:,n), 'ko');
    hold on
    plot(MAST, r_H_re(:,n), 'y-');
    ylim([-0.01 0.01]);
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Tmax (K)','FontSize',10)
    
    subplot(2,3,5);
    n = 5;
    plot(MAST, H_re(:,n), 'ko');
    hold on
    plot(MAST, r_H_re(:,n), 'y-');
    ylim([0 2000]);
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting exponent c','FontSize',10)

    subplot(2,3,6);
    n = 6;
    plot(MAST, H_re(:,n), 'ko');
    hold on
    plot(MAST, r_H_re(:,n), 'y-');
    ylim([0 2000]);
    title('P(6) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter','FontSize',10)

    sgtitle('Huang (Respiration) Parameters vs. Mean Annual Soil Temperature (MAST)');

% MMRT
figure;
    subplot(2,3,1);
    n = 1;
    plot(MAST, M_re(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_M_re(:,n), 'c-');
    ylim([0 14e6]);
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,3,2);
    n = 2;
    plot(MAST, M_re(:,n), 'ko');
    hold on
    plot(MAST, r_M_re(:,n), 'c-');
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,3,3);
    n = 3;
    plot(MAST, M_re(:,n), 'ko');
    hold on
    plot(MAST, r_M_re(:,n), 'c-');
    ylim([-500 8000]);
    title('P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('T0 (K)','FontSize',10)
    
    subplot(2,3,4);
    n = 4;
    plot(MAST, M_re(:,n), 'ko');
    hold on
    plot(MAST, r_M_re(:,n), 'c-');
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)

    subplot(2,3,5);
    n = 5;
    plot(MAST, M_re(:,n), 'ko');
    hold on
    plot(MAST, r_M_re(:,n), 'c-');
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter','FontSize',10)
    
    sgtitle('MMRT (Respiration) Parameters vs. Mean Annual Soil Temperature (MAST)');

% MMRTD
figure;
    subplot(2,4,1);
    n = 1;
    plot(MAST, D_re(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([2e6 8e6]);
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,4,2);
    n = 2;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([-5000 -1000]);
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,4,3);
    n = 3;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([150 1000])
    title('MMRT_D (Growth) P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('T0 (K)','FontSize',10)
    
    subplot(2,4,4);
    n = 4;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([-0.5e4 1.5e4])
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,5);
    n = 5;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Difference in Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,4,6);
    n = 6;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([-100 20]);
    title('P(6) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Difference in Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,4,7);
    n = 7;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([0 80]);
    title('P(7) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter a','FontSize',10)
    
    subplot(2,4,8);
    n = 8;
    plot(MAST, D_re(:,n), 'ko');
    hold on
    plot(MAST, r_D_re(:,n), 'm-');
    ylim([0 150]);
    title('P(8) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter b','FontSize',10)

    sgtitle('MMRT_D (Respiration) Parameters vs. Mean Annual Soil Temperature (MAST)');

% Combined
% MMRT_Dcomb
figure;
    subplot(2,4,1);
    n = 1;
    plot(MAST, D_co(:,n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([3e6 7e6]);
    title('P(1) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,4,2);
    n = 2;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([-6000 -500]);
    title('P(2) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,4,3);
    n = 3;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    title('MMRT_D (Growth) P(3) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('T0 (K)','FontSize',10)
    
    subplot(2,4,4);
    n = 4;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([0 12000]);
    title('P(4) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Heat capacity (kJ/mol/K)','FontSize',10)
    
    subplot(2,4,5);
    n = 5;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([3e4 12e4]);
    title('P(5) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Difference in Enthalpy (kj/mol)','FontSize',10)
    
    subplot(2,4,6);
    n = 6;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([-200 300]);
    title('P(6) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Difference in Entropy (J/mol/K)','FontSize',10)
    
    subplot(2,4,7);
    n = 7;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([0 50]);
    title('P(7) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter a','FontSize',10)
    
    subplot(2,4,8);
    n = 8;
    plot(MAST, D_co(:,n), 'ko');
    hold on
    plot(MAST, r_D_co(:,n), 'g-');
    ylim([50 250]);
    title('P(8) vs MAST');
    xlabel('Temperature (K)','FontSize',10)
    ylabel('Fitting parameter b','FontSize',10)

    sgtitle('MMRT_Dcomb Parameters vs. Mean Annual Soil Temperature (MAST)');

% Table exporting
modelNames = {'Walker','Huang','MMRT','MMRTD', 'MMRTD_comb'};
metricNames = {'R2','RMSE','AIC'};

% Statistics
% Column order: Walker, Huang, MMRT, MMRTD & Row order: R2, RMSE, AIC
summary_gr   = [median(R2_gr(:,1)), median(R2_gr(:,2)), median(R2_gr(:,3)), median(R2_gr(:,4)),median(R2_gr(:,5));
                median(RMSE_gr(:,1)), median(RMSE_gr(:,2)), median(RMSE_gr(:,3)), median(RMSE_gr(:,4)),median(RMSE_gr(:,5));
                median(AIC_gr(:,1)), median(AIC_gr(:,2)), median(AIC_gr(:,3)), median(AIC_gr(:,4)),median(AIC_gr(:,5))];
summary_re   = [median(R2_re(:,1)), median(R2_re(:,2)), median(R2_re(:,3)), median(R2_re(:,4)),median(R2_re(:,5));
                median(RMSE_re(:,1)), median(RMSE_re(:,2)), median(RMSE_re(:,3)), median(RMSE_re(:,4)),median(RMSE_re(:,5));
                median(AIC_re(:,1)), median(AIC_re(:,2)), median(AIC_re(:,3)), median(AIC_re(:,4)),median(AIC_re(:,5))];
summary_co   = [median(R2_co(:,1)),median(R2_co(:,2)),median(R2_co(:,3)),median(R2_co(:,4)),median(R2_co(:,5));
                median(RMSE_co(:,1)),median(RMSE_co(:,2)),median(RMSE_co(:,3)),median(RMSE_co(:,4)),median(RMSE_co(:,5));
                median(AIC_co(:,1)),median(AIC_co(:,2)),median(AIC_co(:,3)),median(AIC_co(:,4)),median(AIC_co(:,5));];  % one-step combined

summary_CUE =  [median(R2_CUE(:,1)), median(R2_CUE(:,2)), median(R2_CUE(:,3)), median(R2_CUE(:,4)),median(R2_CUE(:,5));
                median(RMSE_CUE(:,1)), median(RMSE_CUE(:,2)), median(RMSE_CUE(:,3)), median(RMSE_CUE(:,4)),median(RMSE_CUE(:,5));
                median(AIC_CUE(:,1)), median(AIC_CUE(:,2)), median(AIC_CUE(:,3)), median(AIC_CUE(:,4)),median(AIC_CUE(:,5))];

T_gr = table(summary_gr(1,:)', summary_gr(2,:)', summary_gr(3,:)', ...
             'VariableNames',metricNames, 'RowNames', modelNames);
T_re = table(summary_re(1,:)', summary_re(2,:)', summary_re(3,:)', ...
             'VariableNames',metricNames, 'RowNames', modelNames);
T_co = table(summary_co(1,:)', summary_co(2,:)', summary_co(3,:)', ...
             'VariableNames',{'R2','RMSE','AIC'},'RowNames',{'Walker','Huang','MMRT','MMRTD','MMRTD_one'});
T_CUE = table(summary_CUE(1,:)', summary_CUE(2,:)', summary_CUE(3,:)', ...
             'VariableNames',{'R2','RMSE','AIC'},'RowNames',{'Walker','Huang','MMRT','MMRTD','MMRTD_one'});

writetable(T_gr, 'GrowthMetrics.xlsx', 'WriteRowNames', true);
writetable(T_re, 'RespirationMetrics.xlsx', 'WriteRowNames', true);
writetable(T_co, 'OneStepMetrics.xlsx', 'WriteRowNames', true);
writetable(T_CUE, 'CUEMetrics.xlsx', 'WriteRowNames', true);

% Spearman
colNames = {'Rho','Pval'};

Spearman_Walker_gr = table(Walker_Spearman_gr(:,1), Walker_Spearman_gr(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7'});
Spearman_Huang_gr = table(Huang_Spearman_gr(:,1), Huang_Spearman_gr(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5'});
Spearman_MMRT_gr = table(MMRT_Spearman_gr(:,1), MMRT_Spearman_gr(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4'});
Spearman_MMRTD_gr = table(MMRTD_Spearman_gr(:,1), MMRTD_Spearman_gr(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6'});

Spearman_Walker_re = table(Walker_Spearman_re(:,1), Walker_Spearman_re(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7','P8'});
Spearman_Huang_re = table(Huang_Spearman_re(:,1), Huang_Spearman_re(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6'});
Spearman_MMRT_re = table(MMRT_Spearman_re(:,1), MMRT_Spearman_re(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5'});
Spearman_MMRTD_re = table(MMRTD_Spearman_re(:,1), MMRTD_Spearman_re(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7','P8'});

Spearman_co = table(MMRTD_Spearman_co(:,1), MMRTD_Spearman_co(:,2), ...
    'VariableNames', colNames, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7','P8'});

writetable(Spearman_Walker_gr,'Spearman_Walker_gr.xlsx','WriteRowNames',true);
writetable(Spearman_Huang_gr,'Spearman_Huang_gr.xlsx','WriteRowNames',true);
writetable(Spearman_MMRT_gr,'Spearman_MMRT_gr.xlsx','WriteRowNames',true);
writetable(Spearman_MMRTD_gr,'Spearman_MMRTD_gr.xlsx','WriteRowNames',true);

writetable(Spearman_Walker_re,'Spearman_Walker_re.xlsx','WriteRowNames',true);
writetable(Spearman_Huang_re,'Spearman_Huang_re.xlsx','WriteRowNames',true);
writetable(Spearman_MMRT_re,'Spearman_MMRT_re.xlsx','WriteRowNames',true);
writetable(Spearman_MMRTD_re,'Spearman_MMRTD_re.xlsx','WriteRowNames',true);

writetable(Spearman_co,'Spearman_co.xlsx','WriteRowNames',true);

% Pattern score
colNames2 = {'Pattern score'};

Pat_Walker_gr = table(Walker_Pat_gr(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7'});
Pat_Huang_gr = table(Huang_Pat_gr(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5'});
Pat_MMRT_gr = table(MMRT_Pat_gr(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4'});
Pat_MMRTD_gr = table(MMRTD_Pat_gr(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6'});

Pat_Walker_re = table(Walker_Pat_re(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7','P8'});
Pat_Huang_re = table(Huang_Pat_re(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6'});
Pat_MMRT_re = table(MMRT_Pat_re(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5'});
Pat_MMRTD_re = table(MMRTD_Pat_re(:,1),...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7','P8'});

Pat_co = table(MMRTD_Pat_co(:,1), ...
    'VariableNames', colNames2, ...
    'RowNames', {'P1','P2','P3','P4','P5','P6','P7','P8'});

writetable(Pat_Walker_gr,'Pat_Walker_gr.xlsx','WriteRowNames',true);
writetable(Pat_Huang_gr,'Pat_Huang_gr.xlsx','WriteRowNames',true);
writetable(Pat_MMRT_gr,'Pat_MMRT_gr.xlsx','WriteRowNames',true);
writetable(Pat_MMRTD_gr,'Pat_MMRTD_gr.xlsx','WriteRowNames',true);

writetable(Pat_Walker_re,'Pat_Walker_re.xlsx','WriteRowNames',true);
writetable(Pat_Huang_re,'Pat_Huang_re.xlsx','WriteRowNames',true);
writetable(Pat_MMRT_re,'Pat_MMRT_re.xlsx','WriteRowNames',true);
writetable(Pat_MMRTD_re,'Pat_MMRTD_re.xlsx','WriteRowNames',true);

writetable(Pat_co,'Pat_MMRTD_co.xlsx','WriteRowNames',true);

% MMRT symmetrical function
function [P] = Regression_Tmodel_s1(Tp,Yp,eq,P0) %to call the minimising function
    errfh = @(P,x,z) sum((z(1:7) - (eq(x(1:7),P))).^2)  ...        % main error
                     + 10^6 * sum(abs(imag(eq(x(:),P))));          % penalty for imaginary parts
    Tp_sym = [Tp(1:6) Tp(6)+linspace(5,35,6)'];                    % Making sure the initial guess of MMRT is symmetrical 
    Yp_sym = [Yp(1:6) flip(Yp(1:6))];
    opts = optimset('MaxIter',5000000,'MaxFunEval',5000000,'TolFun',10e-16,'TolX',10e-16,'TolConSQP',10e-16,'Display','off'); % solver options definition
    [P] = fminsearch(errfh,P0,opts,Tp_sym,Yp_sym);%solver function
end

% Huang and Walker function
function [P] = Regression_Tmodel_s2(Tp,Yp,eq,P0)
    errfh = @(P,x,z) sum((z(:) - (eq(x(:),P))).^2) ...  
                     + 10^6 * sum(abs(imag(eq(x(:),P)))); 
    opts = optimset('MaxIter',5000000,'MaxFunEval',5000000,'TolFun',10e-16,'TolX',10e-16,'TolConSQP',10e-16,'Display','off');
    [P] = fminsearch(errfh,P0,opts,Tp,Yp);
end


%MMRTD One step calibration
function [P] = Regression_Tmodel_comb(Tp1,Yp1,eqG,Tp2,Yp2,eqR,P0)
    Tp = [Tp1;Tp2];
    lenG = length(Tp1);
    Yp = [Yp1;Yp2];
    errfh = @(P,x,z) sum([(z(1:lenG)'-eqG(x(1:lenG)',P))./z(1:lenG)' (z(lenG+1:end)'-eqR(x(lenG+1:end)',P))./z(lenG+1:end)'].^2) ...
                      + 0* sum(((P(1:6)-P0(1:6)) ./P0(1:6)) .^2) ...    
                      + 1* sum(((P(7:8)-P0(7:8)) ./P0(7:8)) .^2);            
    opts = optimset('MaxIter',500000,'MaxFunEval',500000,'TolFun',10e-16,'TolX',10e-16,'TolConSQP',10e-16,'Display','off');
    [P] = fmincon(errfh,P0,[],[],[],[],[-inf -inf -inf -inf -inf -inf 10 10],[],[],opts,Tp,Yp);            
end

function [R2,RMSE,AICc,BIC] = funModEval(Yp,Ymod,nPar)
    RSS = sum((Yp-Ymod).^2);            % seperate as it is used in R2 and RMSE            
    R2 = 1 - RSS/sum((Yp-mean(Yp)).^2);             
    nData = length(Yp);                 % number of datapoints
    RMSE = sqrt(RSS/nData); 
    AIC = nData*log(RSS/nData)+2*(nPar+1);
    AICc = AIC + 2*nPar*(nPar+1)/(nData-nPar-1);
    if nData-nPar-1<=0
        pa = -1/(nData-nPar-2);
        AICc = AIC + 2*nPar*(nPar+1)/pa;
    end
    BIC = nData*log(RSS/nData)+nPar*log(nData);
end

function Spearman_mat = spearmanMatrix(MAST, ParamMat)
    nParams = size(ParamMat,2);
    Rho = zeros(nParams,1);
    Pval = zeros(nParams,1);
    
    for p = 1:nParams
        [Rho(p), Pval(p)] = corr(MAST, ParamMat(:,p), 'Type', 'Spearman');
    end
  
    Spearman_mat = [Rho(:), Pval(:)];
end
