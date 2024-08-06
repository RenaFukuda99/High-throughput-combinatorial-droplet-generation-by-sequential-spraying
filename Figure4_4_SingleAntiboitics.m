%% Single Antibiotic 
% Read in single antiboitic droplet data and form dose response curves. 
% Relies on having generated fluorescence calibration curves and 
% growth normalization curve.

% Choose filename from "Reading Images"
filename = "my_filename.xlsx"

% Initialize arrays
h_comb_exp=[];
end_comb_exp=[];
F_conc_comb=[];
R_conc_comb=[];
indices=[];
v_comb_exp = [];

contact_angle = 116;

% Read in data from every sheet/droplet image
for i=1:10 % adjust to number of images/sheets
    
    % Read in sheet and compute height, volume
    T= readtable(filename,'Sheet',i);
    Blue_300_vect = T.BlueInt_300; %If working with solution 2 (linezolid, tobramycin), read in red intensity instead
    h=[T.diameters] ./(2*sind(contact_angle)); 
    v = (1/6)*pi.*(h*3.24).*((h*3.24).^2 + 3*(T.diameters*3.24/2).^2);

    % Initialize arrays 
    R_conc=[];
    F_conc = [];

    % Compute size threshold and background intensity
    mask = T.diameters*3.24 >= 150;
    null_int = 47.3850178607331 + 2.92753870925358*h; %background intensity from positive controls

    % Compute concentration using fluorescent calibration
    for j=1:height(T)
        if Blue_300_vect(j) >= null_int(j) 
            F_conc(j) = (Blue_300_vect(j)) ./(214.0163912*h(j)) ; % value from fluorescent calibration
        else
            F_conc(j) = 0;
        end 
    end

    R_conc(j)=0;
    F_conc(F_conc <= 0) = 0;
    
    % Compute normalized growth from positive controls
    theo_growth =   934.055847255915 + 78.7552567133320*h;
    endpoints = table2array(T(:,17)) - table2array(T(:,6));
    norm_growth=endpoints./theo_growth;
    norm_growth(norm_growth <= 0) = 0;
    
    % Save to vectors
    h_comb_exp=[h_comb_exp; h(mask)];
    v_comb_exp=[v_comb_exp; v(mask)];
    end_comb_exp=[end_comb_exp; norm_growth(mask)];
    F_conc_comb=[F_conc_comb; F_conc(mask)'];
    R_conc_comb=[R_conc_comb; R_conc(mask)'];
end 

%% Dose Response Curve

%Sort data in diameter ascentidng order
[F_conc_comb, a_order] = sort(F_conc_comb);
end_comb_exp = end_comb_exp(a_order,:);

% Fit hill function
fitfun = fittype( @(K,n,x) 1 ./(1+((x./K).^n)),'independent','x', 'dependent', 'y', 'coefficients', {'K','n'} );
options.Lower = [0 0];
options.Upper = [1 8];
[fitted_curve,gof] = fit(F_conc_comb,end_comb_exp,fitfun) %,'StartPoint',x0

% Save the coeffiecient values for a,b,c and d in a vector
coeffvals = coeffvalues(fitted_curve);

% Plot droplet dose response data and fitten Hill curve
figure
scatter(F_conc_comb*500, end_comb_exp/0.7,20,0.000001*v_comb_exp,'MarkerEdgeAlpha',0.7)
hold on
plot(F_conc_comb*500,fitted_curve(F_conc_comb)/0.7,'k')
hold off
title('Random title', 'Color','none')
set(gca,'xscale','log')
ylim([0,1.5])
xlim([10^-1, 10^3])
colormap('parula')
clim([0,20])
xlim([0.1 1000])
