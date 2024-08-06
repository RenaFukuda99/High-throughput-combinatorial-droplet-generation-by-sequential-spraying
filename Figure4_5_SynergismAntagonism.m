%% Synergism/Antagonism Assay
% Read in intensity and droplet information for synergism/antagnosim
% assays. Relies on having generated fluorescence calibration curves, 
% growth normalization curve, and indivudal dose-response curves.

% Read in relevant file
filename = "my_filename.xlsx";

% Initialize arrays
h_comb_exp=[];
end_comb_exp=[];
F_conc_comb=[];
R_conc_comb=[];
indices=[];
contact_angle = 116;

% Read in droplets from each sheet
for i=1:10 % change upper bound to number of sheets 
    
    % Read in data
    T= readtable(filename,'Sheet',i);
    h=[T.diameters] ./(2*sind(contact_angle)); 
    Blue_300_vect = T.BlueInt_300;
    Red_300_vect = T.RedInt_300;

    % Initialize arrays
    R_conc=[];
    F_conc = [];

    % Set mask for size-based filtering
    mask = T.diameters*3.24 >= 150; %3.24 is pixel to um conversion
    
    % Calculate fluorescent concentrations using values from background
    % values and fluorescent calibration curves
    for j=1:height(T)
        if Blue_300_vect(j) >= 5
             F_conc(j) = (Blue_300_vect(j)-40.8294975215378) ./(h(j)*527.589687137501) ;
             R_conc(j) = (Red_300_vect(j)) ./(252.986370219884*h(j)) ;
        else
            F_conc(j) = 0;
             R_conc(j) = (Red_300_vect(j)) ./(252.986370219884*h(j)) ;
        end 
    end
    
    F_conc(F_conc <= 0) = 0;
    R_conc(R_conc <= 0) = 0;
    
    % Calculate theoretical growth from positive control curve
    theo_growth = (0 + 58.439946980099030*h);
    endpoints = table2array(T(:,13)) - table2array(T(:,6));
    norm_growth=endpoints./theo_growth;
    norm_growth(norm_growth <= 0) = 0;
    
    % Save data to vectors
    h_comb_exp=[h_comb_exp; h(mask)];
    end_comb_exp=[end_comb_exp; norm_growth(mask)];
    F_conc_comb=[F_conc_comb; F_conc(mask)'];
    R_conc_comb=[R_conc_comb; R_conc(mask)'];
end 

%% Plot pairwise heatmap
mask_conc = F_conc_comb <= 1 & R_conc_comb <= 1;

% Discretize data (set bins), since data is being shown in a heatmap 
edges_F = linspace(0,1,11);
values_F = edges_F(2:end);
edges_R = linspace(0,1,11);
values_R = edges_R(2:end);
Vanco_conc = discretize(F_conc_comb(mask_conc),edges_F, values_F);
Tobra_conc = discretize(R_conc_comb(mask_conc), edges_R, values_R);

% Create bin labels
EFreal=char.empty;
ERreal=char.empty;
for i=1:10
    EFreal{i} = strcat('[', num2str(edges_F(i)), ', ', num2str(edges_F(i+1)), ']');
    ERreal{i} = strcat('[', num2str(edges_R(i)), ', ', num2str(edges_R(i+1)), ']');
end

% Plot heatmap of concentrations against normalized intensities
end_comb_exp = end_comb_exp(mask_conc);
tbl = table(Vanco_conc,Tobra_conc,end_comb_exp);
hm=heatmap(tbl,'Vanco_conc','Tobra_conc','ColorVariable','end_comb_exp');
hm.YDisplayData = flipud(hm.YDisplayData);
hm.clim([0 1])

% Compute number of droplets per bin
disc_F = discretize(F_conc_comb(mask_conc),linspace(min(F_conc_comb(mask_conc)),1,11));
disc_R = discretize(R_conc_comb(mask_conc),linspace(min(R_conc_comb(mask_conc)),1,11));

% Plot droplets per bin as heatmap 
tbl2=table(disc_F,disc_R);
hdisc=heatmap(tbl2,'disc_F','disc_R')

% Bliss score calculation 
yf = 1 - (1 ./(1+((tbl.Vanco_conc./ 0.5276).^2.781))); % K, n values are from individual dose response curves from single antibiotic controls
yr = 1 - (1 ./(1+((tbl.Tobra_conc./0.99).^6))); % K, n values are from individual dose response curves from single antibiotic controls
perc_inhib = 1-end_comb_exp;
bliss_rate = (yr) + (yf) - (yr).*(yf);
bls_score = perc_inhib - bliss_rate;

% Plot bliss score as heatmap 
tbl_bls = table(Vanco_conc,Tobra_conc,bls_score);
h_bls=heatmap(tbl_bls,'Vanco_conc','Tobra_conc','ColorVariable','bls_score');
h_bls.YDisplayData = flipud(h_bls.YDisplayData);
h_bls.clim([-1 1])
colormap cool

% Save data in table
my_table =table(end_comb_exp, perc_inhib, yr,yf, bliss_rate, bls_score);
