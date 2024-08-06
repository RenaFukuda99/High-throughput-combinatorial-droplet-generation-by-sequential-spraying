%% Fluorescent Calibration
% Calibrate fluorophore intensity to diameter relationship 

% Retrieve relevant file (droplets with only spray 2 or spray 3) 
filename = "my_filename.xlsx";

%% Fluorescein
% Initialize vectors
B_h_comb_F = [];
B_B30_int_comb = [];
B_B300_int_comb = [];
B_R30_int_comb = [];
B_R300_int_comb = [];

% Retrieve data from fluorescein droplets
for i=1:10
    contact_angle = 116;
    T= readtable(filename,'Sheet',i);
    h_now=[T.diameters] ./(2*sind(contact_angle));
    B_B30_int_comb = [B_B30_int_comb; T.BlueInt_30];
    B_B300_int_comb = [B_B300_int_comb; T.BlueInt_300];
    B_R30_int_comb = [B_R30_int_comb; T.RedInt_30];
    B_R300_int_comb = [B_R300_int_comb; T.RedInt_300];
    B_h_comb_F = [B_h_comb_F; h_now];
end 

% Plot intensity vs diameter
figure
tcl = tiledlayout(2,2);
%Plot Blue 30
nexttile
scatter(B_h_comb_F, B_B30_int_comb);
dlm = fitlm(B_h_comb_F,B_B30_int_comb, intercept=false);
hold on
plot(dlm);
B_slope_B30 = dlm.Coefficients.Estimate;
title('Blue 30')

%Plot Blue 300
nexttile
scatter(B_h_comb_F, B_B300_int_comb);
dlm = fitlm(B_h_comb_F(B_h_comb_F<=8000),B_B300_int_comb(B_h_comb_F<=8000), intercept=false);
hold on
plot(dlm);
B_slope_B300 = dlm.Coefficients.Estimate;
title('Blue 300')

%Plot Red 30
nexttile
scatter(B_h_comb_F, B_R30_int_comb);
dlm = fitlm(B_h_comb_F,B_R30_int_comb, intercept=false);
hold on
plot(dlm);
B_slope_R30 = dlm.Coefficients.Estimate;
title('Red 30')

% Plot Red 300
nexttile
scatter(B_h_comb_F, B_R300_int_comb);
dlm = fitlm(B_h_comb_F,B_R300_int_comb, intercept=false);
hold on
plot(dlm);
B_slope_R300 = dlm.Coefficients.Estimate;
title('Red 300')

title(tcl,'Fluorescein Calibration')

%% Alexa Fluor
% Initialize vectors
R_h_comb_F = [];
R_B30_int_comb = [];
R_B300_int_comb = [];
R_R30_int_comb = [];
R_R300_int_comb = [];

% Retrieve data from fluorescein droplets
for i=11:20
    contact_angle = 116;
    T= readtable(filename,'Sheet',i);
    h_now=[T.diameters] ./(2*sind(contact_angle));
    R_B30_int_comb = [R_B30_int_comb; T.BlueInt_30];
    R_B300_int_comb = [R_B300_int_comb; T.BlueInt_300];
    R_R30_int_comb = [R_R30_int_comb; T.RedInt_30];
    R_R300_int_comb = [R_R300_int_comb; T.RedInt_300];
    R_h_comb_F = [R_h_comb_F; h_now];
end 

% Plot intensity vs diameter
figure
%Plot Blue 30
tcl2 = tiledlayout(2,2);

nexttile
scatter(R_h_comb_F, R_B30_int_comb);
dlm = fitlm(R_h_comb_F(R_B30_int_comb<=400),R_B30_int_comb(R_B30_int_comb<=400), intercept=false);
hold on
plot(dlm);
R_slope_B30 = dlm.Coefficients.Estimate;
title('Blue 30')

%Plot Blue 300
nexttile
scatter(R_h_comb_F, R_B300_int_comb);
dlm = fitlm(R_h_comb_F(R_B300_int_comb<=8000),R_B300_int_comb(R_B300_int_comb<=8000), intercept=false);
hold on
plot(dlm);
R_slope_B300 = dlm.Coefficients.Estimate;
title('Blue 300')

%Plot Red 30
nexttile
scatter(R_h_comb_F, R_R30_int_comb);
dlm = fitlm(R_h_comb_F,R_R30_int_comb, intercept=false);
hold on
plot(dlm);
R_slope_R30 = dlm.Coefficients.Estimate;
title('Red 30')

% Plot Red 300
nexttile
scatter(R_h_comb_F, R_R300_int_comb);
dlm = fitlm(R_h_comb_F,R_R300_int_comb, intercept=false);
hold on
plot(dlm);
R_slope_R300 = dlm.Coefficients.Estimate;
title('Red 300')

title(tcl2,'Alexa Fluor Calibration')

% Export coefficients to file
T=table(B_slope_B30, B_slope_B300, B_slope_R30, B_slope_R300, R_slope_B30, R_slope_B300, R_slope_R30, R_slope_R300);
writetable(T,'Coeff.xlsx','Sheet',1);
