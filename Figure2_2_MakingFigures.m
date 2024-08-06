%% Coefficient of Variation 
% Plotting coefficient of variation and percent deviation for each
% fluorophore. Requires reading in images and exporting out tables. 
clc;clear;

%% Read in data and assign concentration
filename_Af = "AlexaFluor.xlsx"; 
filename_Fluor = "Fluorescein.xlsx";

% Alexa fluor - import different concentrations and intensities
AF_D1 = readtable(filename_Af, sheet = "D1");
AF_D1.Conc = repelem(1,height(AF_D1))';
AF_D2 = readtable(filename_Af, sheet = "D2");
AF_D2.Conc = repelem(0.5,height(AF_D2))';
AF_D3 = readtable(filename_Af, sheet = "D3");
AF_D3.Conc = repelem(0.25,height(AF_D3))';
AF_D4 = readtable(filename_Af, sheet = "D4");
AF_D4.Conc = repelem(0.125,height(AF_D4))';
AF_D5 = readtable(filename_Af, sheet = "D5");
AF_D5.Conc = repelem(0.0625,height(AF_D5))';
AF_D6 = readtable(filename_Af, sheet = "D6");
AF_D6.Conc = repelem(0.0312,height(AF_D6))';
AF_D7 = readtable(filename_Af, sheet = "D7");
AF_D7.Conc = repelem(0.0156,height(AF_D7))';

% Fluorescein - import different concentrations and intensities
F_D1 = readtable(filename_Fluor, sheet = "D1");
F_D1.Conc = repelem(1,height(F_D1))';
F_D2 = readtable(filename_Fluor, sheet = "D2");
F_D2.Conc = repelem(0.5,height(F_D2))';
F_D3 = readtable(filename_Fluor, sheet = "D3");
F_D3.Conc = repelem(0.25,height(F_D3))';
F_D4 = readtable(filename_Fluor, sheet = "D4");
F_D4.Conc = repelem(0.125,height(F_D4))';
F_D5 = readtable(filename_Fluor, sheet = "D5");
F_D5.Conc = repelem(0.0625,height(F_D5))';
F_D6 = readtable(filename_Fluor, sheet = "D6");
F_D6.Conc = repelem(0.0312,height(F_D6))';
F_D7 = readtable(filename_Fluor, sheet = "D7");
F_D7.Conc = repelem(0.0156,height(F_D7))';

F_total = [F_D1; F_D2; F_D3; F_D4; F_D5; F_D6; F_D7];
AF_total = [AF_D1; AF_D2; AF_D3; AF_D4; AF_D5; AF_D6; AF_D7; AF_D8 ];

%clearvars -except F_total AF_total

%% Alexa Fluor - Scatter plots of intensity against diameter with linear fits
% Diameters have 3.24 calibration from pixels to um
% Fit linear curvve for each concentration 

scatter(AF_D1.diameters*3.24, AF_D1.Int_3,'o',"MarkerEdgeColor",[255/255 74/255 74/255]);
dlm1 = fitlm(AF_D1{AF_D1.diameters >50, 'diameters'}*3.24,AF_D1{AF_D1.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(1) = dlm1.Coefficients.Estimate;

scatter(AF_D2.diameters*3.24, AF_D2.Int_3,'o',"MarkerEdgeColor",[255/255 140/255 74/255]);
dlm1 = fitlm(AF_D2{AF_D2.diameters >50, 'diameters'}*3.24,AF_D2{AF_D2.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(2) = dlm1.Coefficients.Estimate;

scatter(AF_D3.diameters*3.24, AF_D3.Int_3,'o',"MarkerEdgeColor",[255/255 207/255 74/255]);
dlm1 = fitlm(AF_D3{AF_D3.diameters >50, 'diameters'}*3.24,AF_D3{AF_D3.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(3) = dlm1.Coefficients.Estimate;

scatter(AF_D4.diameters*3.24, AF_D4.Int_3,'o',"MarkerEdgeColor",[201/255 255/255 74/255]);
dlm1 = fitlm(AF_D4{AF_D4.diameters >50, 'diameters'}*3.24,AF_D4{AF_D4.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(4) = dlm1.Coefficients.Estimate;

scatter(AF_D5.diameters*3.24, AF_D5.Int_3,'o',"MarkerEdgeColor",[95/255 222/255 67/255]);
dlm1 = fitlm(AF_D5{AF_D5.diameters >50, 'diameters'}*3.24,AF_D5{AF_D5.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(5) = dlm1.Coefficients.Estimate;

scatter(AF_D6.diameters*3.24, AF_D6.Int_3,'o',"MarkerEdgeColor",[67/255 222/255 196/255]);
dlm1 = fitlm(AF_D6{AF_D6.diameters >50, 'diameters'}*3.24, AF_D6{AF_D6.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(6) = dlm1.Coefficients.Estimate;

scatter(AF_D7.diameters*3.24, AF_D7.Int_3,'o',"MarkerEdgeColor",[155/255 67/255 222/255]);
dlm1 = fitlm(AF_D7{AF_D7.diameters >50, 'diameters'}*3.24,AF_D7{AF_D7.diameters >50, 'Int_3'}, intercept=false);
hold on
plot(dlm1, 'Marker', 'none');
slope_AF(7) = dlm1.Coefficients.Estimate;

xlabel('Diameter (um)')
ylabel('Intensity (A.U.)')
title("Alexa Fluor")

% Scatterplot of slopes against diameter (show proportionality of
% concentration against intensity)
figure;
scatter([1 0.5 0.25 0.125 0.0625 0.0312 0.0156], slope_AF,[], 'filled');
hold on
dlm_slope_AF = fitlm([1 0.5 0.25 0.125 0.0625 0.0312 0.0156], slope_AF);
plot(dlm_slope_AF, 'Marker', 'none')
ylim([0 6])

%% Alexa Fluor - Bar plots of CoV
% Calculate coefficient of variation for each concentration
diam_thresh = 50;

proj_c1 = AF_D1{AF_D1.diameters >diam_thresh, 'Int_3'} ./ AF_D1{AF_D1.diameters >diam_thresh, 'diameters'};
cov_AF(1) = std(proj_c1)/mean(proj_c1);

proj_c2 = AF_D2{AF_D2.diameters >diam_thresh, 'Int_3'} ./ AF_D2{AF_D2.diameters >diam_thresh, 'diameters'};
cov_AF(2) = std(proj_c2)/mean(proj_c2);

proj_c3 = AF_D3{AF_D3.diameters >diam_thresh, 'Int_3'} ./ AF_D3{AF_D3.diameters >diam_thresh, 'diameters'};
cov_AF(3) = std(proj_c3)/mean(proj_c3);

proj_c4 = AF_D4{AF_D4.diameters >diam_thresh, 'Int_3'} ./ AF_D4{AF_D4.diameters >diam_thresh, 'diameters'};
cov_AF(4) = std(proj_c4)/mean(proj_c4);

proj_c5 = AF_D5{AF_D5.diameters >diam_thresh, 'Int_3'} ./ AF_D5{AF_D5.diameters >diam_thresh, 'diameters'};
cov_AF(5) = std(proj_c5)/mean(proj_c5);

proj_c6 = AF_D6{AF_D6.diameters >diam_thresh, 'Int_3'} ./ AF_D6{AF_D6.diameters >diam_thresh, 'diameters'};
cov_AF(6) = std(proj_c6)/mean(proj_c6);

proj_c7 = AF_D7{AF_D7.diameters >diam_thresh, 'Int_3'} ./ AF_D7{AF_D7.diameters >diam_thresh, 'diameters'};
cov_AF(7) = std(proj_c7)/mean(proj_c7);

proj_c8 = AF_D8{AF_D8.diameters >diam_thresh, 'Int_3'} ./ AF_D8{AF_D8.diameters >diam_thresh, 'diameters'};
cov_AF(8) = std(proj_c8)/mean(proj_c8);

mean(cov_AF(1:8))

% Bar plot
figure;
b = bar(cov_AF(1:8)*100);
hold on
title("Alexa Fluor")
ylabel("Coefficient of Variation (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0 20]);
b.EdgeColor = 'none';

%% Fluorescein - Scatter plots with dlm  

figure;
scatter(F_D1.diameters*3.24, F_D1.Int_3,'o',"MarkerEdgeColor",[255/255 74/255 74/255]);
dlm1 = fitlm(F_D1{F_D1.diameters >50, 'diameters'}*3.24, F_D1{F_D1.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(1) = dlm1.Coefficients.Estimate(2);

scatter(F_D2.diameters*3.24, F_D2.Int_3,'o',"MarkerEdgeColor",[255/255 140/255 74/255]);
dlm1 = fitlm(F_D2{F_D2.diameters >50, 'diameters'}*3.24,F_D2{F_D2.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(2) = dlm1.Coefficients.Estimate(2);

scatter(F_D3.diameters*3.24, F_D3.Int_3,'o',"MarkerEdgeColor",[255/255 207/255 74/255]);
dlm1 = fitlm(F_D3{F_D3.diameters >50, 'diameters'}*3.24,F_D3{F_D3.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(3) = dlm1.Coefficients.Estimate(2);

scatter(F_D4.diameters*3.24, F_D4.Int_3,'o',"MarkerEdgeColor",[201/255 255/255 74/255]);
dlm1 = fitlm(F_D4{F_D4.diameters >50, 'diameters'}*3.24,F_D4{F_D4.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(4) = dlm1.Coefficients.Estimate(2);

scatter(F_D5.diameters*3.24, F_D5.Int_3,'o',"MarkerEdgeColor",[95/255 222/255 67/255]);
dlm1 = fitlm(F_D5{F_D5.diameters >50, 'diameters'}*3.24,F_D5{F_D5.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(5) = dlm1.Coefficients.Estimate(2);

scatter(F_D6.diameters*3.24, F_D6.Int_3,'o',"MarkerEdgeColor",[67/255 222/255 196/255]);
dlm1 = fitlm(F_D6{F_D6.diameters >50, 'diameters'}*3.24,F_D6{F_D6.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(6) = dlm1.Coefficients.Estimate(2);

scatter(F_D7.diameters*3.24, F_D7.Int_3,'o',"MarkerEdgeColor",[155/255 67/255 222/255]);
dlm1 = fitlm(F_D7{F_D7.diameters >50, 'diameters'}*3.24,F_D7{F_D7.diameters >50, 'Int_3'});
hold on
plot(dlm1, 'Marker', 'none');
slope_F(7) = dlm1.Coefficients.Estimate(2);

xlabel('Diameter (um)')
ylabel('Intensity (A.U.)')
title("Fluorescein")

% Scatter of slopes
figure;
scatter([1 0.5 0.25 0.125 0.0625 0.0312 0.0156], slope_F,[], 'filled');
hold on
dlm_slope_F = fitlm([1 0.5 0.25 0.125 0.0625 0.0312 0.0156], slope_F);
plot(dlm_slope_F, 'Marker', 'none')
ylim([0 10])

%% Fluorescein - Bar plots of CoV
proj_c1 = F_D1{F_D1.diameters >50, 'Int_3'} ./ F_D1{F_D1.diameters >50, 'diameters'};
cov_F(1) = std(proj_c1)/mean(proj_c1);

proj_c2 = F_D2{F_D2.diameters >50, 'Int_3'} ./ F_D2{F_D2.diameters >50, 'diameters'};
cov_F(2) = std(proj_c2)/mean(proj_c2);

proj_c3 = F_D3{F_D3.diameters >50, 'Int_3'} ./ F_D3{F_D3.diameters >50, 'diameters'};
cov_F(3) = std(proj_c3)/mean(proj_c3);

proj_c4 = F_D4{F_D4.diameters >50, 'Int_3'} ./ F_D4{F_D4.diameters >50, 'diameters'};
cov_F(4) = std(proj_c4)/mean(proj_c4);

proj_c5 = F_D5{F_D5.diameters >50, 'Int_3'} ./ F_D5{F_D5.diameters >50, 'diameters'};
cov_F(5) = std(proj_c5)/mean(proj_c5);

proj_c6 = F_D6{F_D6.diameters >50, 'Int_3'} ./ F_D6{F_D6.diameters >50, 'diameters'};
cov_F(6) = std(proj_c6)/mean(proj_c6);

proj_c7 = F_D7{F_D7.diameters >50, 'Int_3'} ./ F_D7{F_D7.diameters >50, 'diameters'};
cov_F(7) = std(proj_c7)/mean(proj_c7);

mean(cov_F)

% Bar plot
figure;
b = bar(cov_F*100);
hold on
title("Fluorescein")
ylabel("Coefficient of Variation (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0 20]);
b.EdgeColor = 'none';

%% Histograms of density and volume
% Compute droplet volume for each droplet
diams = [AF_total.diameters; F_total.diameters]*3.24;
h = diams ./(2*sind(116)); 
v = (1/6)*pi.*(h).*((h).^2 + 3*(diams/2).^2);

d_test = 50;
h_test = d_test/(2*sind(116));
v_test = (1/6)*pi.*(h_test*3.24).*((h_test*3.24).^2 + 3*(d_test*3.24/2).^2);

% Plot histogram of droplet volumes
figure;
[~,edges] = histcounts(log10(v(v >= 0)*0.000001), 25);
edges = [-5:0.2:3]
histogram(v(v >= 0)*0.000001,10.^edges,'FaceColor',[0 0.6 0.5])
set(gca, 'xscale','log')
set(gca, 'yscale','log')
xlabel('Volume (nL)')
ylabel('Counts')
ylim([0 6000])
hold on
xline(v_test*0.000001, 'LineWidth',3)

mean(v((diams/3.24) >= 50))*0.000001

% Compute density
F_total = [F_D1; F_D2; F_D3; F_D4; F_D5; F_D6; F_D7];
AF_total = [AF_D1; AF_D2; AF_D3; AF_D4; AF_D5; AF_D6; AF_D7; AF_D8 ];

F_total = F_total(F_total.diameters >= 50);
AF_total = F_total(F_total.diameters >= 50);

[~,~,ix] = unique(AF_total.Var1);
AF_density = accumarray(ix,1).'

[~,~,ix] = unique(F_total.Var1);
F_density = accumarray(ix,1).'

density = [AF_density F_density];

% Plot density histogram
figure;
histogram(density/(2044*2048*3.24*3.24*(10^-8)),[250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000],'FaceColor',[0 0.6 0.5],'FaceAlpha',0.7)
xlabel('Droplets per cm^2')
ylabel('Counts')

mean(density)/(2044*2048*3.24*3.24*(10^-8))


%% Variation from fit - compute fractional/percent deviation from fit
% Alexa Fluor
perc_dev_D1 = abs((AF_D1.diameters*3.24*slope_AF(1) - AF_D1.Int_3)./(AF_D1.diameters*3.24*slope_AF(1)));
perc_dev_D2 = abs((AF_D2.diameters*3.24*slope_AF(2) - AF_D2.Int_3)./(AF_D2.diameters*3.24*slope_AF(2)));
perc_dev_D3 = abs((AF_D3.diameters*3.24*slope_AF(3) - AF_D3.Int_3)./(AF_D3.diameters*3.24*slope_AF(3)));
perc_dev_D4 = abs((AF_D4.diameters*3.24*slope_AF(4) - AF_D4.Int_3)./(AF_D4.diameters*3.24*slope_AF(4)));
perc_dev_D5 = abs((AF_D5.diameters*3.24*slope_AF(5) - AF_D5.Int_3)./(AF_D5.diameters*3.24*slope_AF(5)));
perc_dev_D6 = abs((AF_D6.diameters*3.24*slope_AF(6) - AF_D6.Int_3)./(AF_D6.diameters*3.24*slope_AF(6)));
perc_dev_D7 = abs((AF_D7.diameters*3.24*slope_AF(7) - AF_D7.Int_3)./(AF_D7.diameters*3.24*slope_AF(7)));

perc_dev_AF = [perc_dev_D1; perc_dev_D2; perc_dev_D3; perc_dev_D4; perc_dev_D5; perc_dev_D6; perc_dev_D7];
AF_noD8 = AF_total(AF_total.Conc ~= 0.0078, :);

figure;
scatter(AF_noD8.diameters(AF_noD8.diameters >= 0)*3.24, perc_dev_AF(AF_noD8.diameters >= 0), ...
    [], "filled", ...
    'MarkerEdgeColor',"none", ...
    "MarkerFaceAlpha", 0.05)
set(gca, 'yscale','log')
hold on

[h_now,I] = sort(AF_noD8.diameters(AF_noD8.diameters >= 0));
unsort = perc_dev_AF(AF_noD8.diameters >= 0);
simple = movavg(unsort(I),'simple',20);
plot(h_now*3.24, simple, "Color", 'k')

% Fluorescein
perc_dev_D1 = abs((F_D1.diameters*3.24*slope_F(1) - F_D1.Int_3)./(F_D1.diameters*3.24*slope_F(1)));
perc_dev_D2 = abs((F_D2.diameters*3.24*slope_F(2) - F_D2.Int_3)./(F_D2.diameters*3.24*slope_F(2)));
perc_dev_D3 = abs((F_D3.diameters*3.24*slope_F(3) - F_D3.Int_3)./(F_D3.diameters*3.24*slope_F(3)));
perc_dev_D4 = abs((F_D4.diameters*3.24*slope_F(4) - F_D4.Int_3)./(F_D4.diameters*3.24*slope_F(4)));
perc_dev_D5 = abs((F_D5.diameters*3.24*slope_F(5) - F_D5.Int_3)./(F_D5.diameters*3.24*slope_F(5)));
perc_dev_D6 = abs((F_D6.diameters*3.24*slope_F(6) - F_D6.Int_3)./(F_D6.diameters*3.24*slope_F(6)));
perc_dev_D7 = abs((F_D7.diameters*3.24*slope_F(7) - F_D7.Int_3)./(F_D7.diameters*3.24*slope_F(7)));

perc_dev_F = [perc_dev_D1; perc_dev_D2; perc_dev_D3; perc_dev_D4; perc_dev_D5; perc_dev_D6; perc_dev_D7];

figure;
scatter(F_total.diameters(F_total.diameters >= 0)*3.24, perc_dev_F(F_total.diameters >= 0), ...
    [], "filled", ...
    'MarkerEdgeColor',"none", ...
    "MarkerFaceAlpha", 0.05)
set(gca, 'yscale','log')
hold on

[h_now,I] = sort(F_total.diameters(F_total.diameters >= 0));
unsort = perc_dev_F(F_total.diameters >= 0);
simple = movavg(unsort(I),'simple',20);
plot(h_now*3.24, simple, "Color", 'k')
