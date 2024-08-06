%% Positive control calibration 
% Used on positive control droplets (no antibiotic or fluorophores)
% Used to obtain normalized growth and background fluorescence values

% Get relevant file (from Figure 4_1_ReadingImages for positive controls)
filename = "my_filenams.xlsx";

% Initialize vectors
h_comb=[];
end_comb=[];
norm_end_comb = [];
contact_angle = 116;

% Loop through sheets to obtain normalized intensity information
for i=1:10
    % Read in table
    T= readtable(filename,'Sheet',i);
    h=[T.diameters] ./(2*sind(contact_angle)); 
    
    % Filter by size
    mask = T.diameters*3.24 >= 150; %3.24 is conversion from pixels to um
    
    % Get endpoint change in intensity
    endpoints = table2array(T(:,6:17)) - repmat(table2array(T(:,6)),1,12);
    
    % Save values to vectors
    h_comb=[h_comb; h];
    end_comb=[end_comb; endpoints];
    norm_end_comb = [norm_end_comb; endpoints(mask,:)./h(mask)];

    clearvars -except slope_F filename h_comb end_comb norm_end_comb contact_angle
end 

% Plot normalized intensity against height
figure
scatter(h_comb(mask), end_comb(mask), 'filled')
hold on
dlm = fitlm(h_comb(mask), end_comb(mask))
plot(dlm);
pos_fit=dlm.Coefficients.Estimate % Record these for later-on scripts
hold on
xlabel('Height of droplet')
ylabel('Increase in OD')
title('Positive Controls')

%% Background droplet fluorescence values
% Retrieve backgroun fluorescent curves

% Initialize vectors
h_comb_F = [];
B30_int_comb = [];
B300_int_comb = [];
R30_int_comb = [];
R300_int_comb = [];

% Read in fluorescent data for positive control droplets
for i=1:10
    contact_angle = 116;
    T= readtable(filename,'Sheet',i);
    h_now=[T.diameters] ./(2*sind(contact_angle));
    B30_int_comb = [B30_int_comb; T.BlueInt_30];
    B300_int_comb = [B300_int_comb; T.BlueInt_300];
    R30_int_comb = [R30_int_comb; T.RedInt_30];
    R300_int_comb = [R300_int_comb; T.RedInt_300];
    h_comb_F = [h_comb_F; h_now];
end 

% Plot fluorescence vs diameter  
figure
tcl = tiledlayout(2,2);
%Plot Blue 30
nexttile
scatter(h_comb_F, B30_int_comb);
dlm = fitlm(h_comb_F,B30_int_comb);
hold on
plot(dlm);
B_slope_B30 = dlm.Coefficients.Estimate;
title('Blue 30')

%Plot Blue 300
nexttile
scatter(h_comb_F, B300_int_comb);
dlm = fitlm(h_comb_F,B300_int_comb);
hold on
plot(dlm);
B_slope_B300 = dlm.Coefficients.Estimate;
title('Blue 300')

%Plot Red 30
nexttile
scatter(h_comb_F, R30_int_comb);
dlm = fitlm(h_comb_F,R30_int_comb);
hold on
plot(dlm);
B_slope_R30 = dlm.Coefficients.Estimate;
title('Red 30')

% Plot Red 300
nexttile
scatter(h_comb_F, R300_int_comb);
dlm = fitlm(h_comb_F,R300_int_comb);
hold on
plot(dlm);
B_slope_R300 = dlm.Coefficients.Estimate;
title('Red 300')

title(tcl,'Blue Channel')

% Export to table
T=table(B_slope_B30, B_slope_B300, B_slope_R30, B_slope_R300);
