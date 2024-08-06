%% Positive control (Figure 3)

clc;clear;

% Get directories
D=dir();
Dark = D(contains({D.name},'4X_Dark.tif'));

% Set up darkfield intensity threshold
thresh_dia = 15000;

% Set up save file
filename = 'FigureThree.xlsx';
positions = [41:2:59]; % XY positions of images

%% Read timelapse darkfield intensities
% If excel file is already generated, can skip this section

for j=1:10
    
    % Obtain structure of relevant images (XY position) 
    Dark_now =  Dark(contains({Dark.name}, strcat('XY',num2str(positions(j),'%02.f'))));

    % Create preliminary mask
    im_dark = imread(strcat(Dark_now(1).folder, '\', Dark_now(1).name));
    im_dark = imresize(im_dark, 2);
    mask = imfill(im_dark >= thresh_dia, 'holes');
    mask = bwareaopen(mask, 10);
    mask_raw = mask;

    % Remove edges to get rid of partial droplets
    mask(1:20,:) = 0; 
    mask(:,1:20) = 0;
    mask(2020:2044,:) = 0;
    mask(:,2020:2048) = 0;

    % Filter preliminary mask
    droplets = regionprops(mask, 'Centroid', 'Area', 'Circularity','MinorAxisLength', 'MajorAxisLength');
    droplets = droplets([droplets.Circularity] > 0.99);
    droplets = droplets([droplets.Area] <= 100000);
    droplets = droplets([droplets.Area] >= 0);

    % Read droplet radii and centers from mask
    diameters_1=([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2;
    centers = reshape([droplets.Centroid], 2, width([droplets.Centroid])/2)';
    radii = diameters_1/2;

%     % Visualize circles
%     imshow(im_dark);
%     hold on
%     viscircles(centers,radii);

    % Make circular mask
    [X Y] = size(im_dark);
    [columnsInImage rowsInImage] = meshgrid(1:Y, 1:X);
    mask_temp = zeros(size(im_dark));

    for k=1:length(droplets)
        centroid_x(k) = droplets(k).Centroid(2);
        centroid_y(k) = droplets(k).Centroid(1);
        circlePixels = (rowsInImage - centroid_x(k)).^2 + (columnsInImage - centroid_y(k)).^2 <= radii(k).^2;
        mask_temp(circlePixels ~= 0) = true;
        mask = mask_temp;
    end

    % Read droplet intensities
    labeled_mask = bwlabel(mask, 8);
    droplets = regionprops(labeled_mask, 'Area', 'Circularity','MinorAxisLength', 'MajorAxisLength','Centroid');
    diameters=([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2;

    % Make smaller circular mask
    [X Y] = size(im_dark);
    [columnsInImage rowsInImage] = meshgrid(1:Y, 1:X);
    radii_sm = 0.5*((([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2)/2);
    filler_mask = zeros(X,Y);

    for k=1:length(droplets)
        centroid_x(k) = droplets(k).Centroid(2);
        centroid_y(k) = droplets(k).Centroid(1);
        circlePixels = (rowsInImage - centroid_x(k)).^2 + (columnsInImage - centroid_y(k)).^2 <= (radii_sm(k)).^2;
        filler_mask(circlePixels ~= 0) = true;
    end
    
    % Label smaller mask
    lab_mask_small = filler_mask.*labeled_mask;

    % Get darkfield intensity information
    dark_int=[];
    
    for i=1:8
        
        % Read darkfield image
        im = imread(strcat(Dark_now(i).folder, '\', Dark_now(i).name));
        im = imresize(im, 2);

        % Read intensity
        im(im >= 40000) = NaN; % Remove droplet edges 
        im = imtranslate(im, [(-1.9 - 1.5*i) 0.3*i]);
        drop_dark = regionprops(lab_mask_small, im, 'MeanIntensity');
        dark_int(:,i)= [drop_dark.MeanIntensity]';
    end 
    
    % Write to table
    T=table(diameters, dark_int);
    writetable(T,filename,'Sheet',(j));

    clearvars -except Dark Dia thresh_dia filename positions
    close all
end

%% Read in data from table
% Initialize arrays
h_comb=[];
end_comb=[];
norm_end_comb=[];
collapsed = [];
contact_angle = 116;

for i=1:10
    % Read table and calculate height
    T= readtable(filename,'Sheet',i);
    dark_array = table2array(T(:,2:9)); 
    h=[T.diameters] ./(2*sind(contact_angle)); 

    % Subtract off minimum intensities
    min_array = min(dark_array,[],2);
    sub_array = dark_array - repmat(min_array,[1 8]); 
    endpoints = sub_array;
    
    % Filter droplets by size, only choose droplets with growth   
    mask_dia = (dark_array(:,8) >= dark_array(:,1)) & (dark_array(:,3) >= dark_array(:,2))  & (dark_array(:,6) >= dark_array(:,5));
    mask = mask_dia & T.diameters*3.24 >= 150 & sub_array(:,1)==0; % 3.24 is conversion of pixels to um

    % Save height, intensity, and normalized intensities to vectors
    h_comb=[h_comb; h(mask)];
    end_comb=[end_comb; endpoints(mask,:)];
    norm_end_comb = [norm_end_comb; endpoints(mask,:)./h(mask)];

    clearvars -except slope_F filename h_comb end_comb collapsed norm_end_comb contact_angle
end 

%% Figures

% Plotting raw intensity for each droplet
figure
subplot(1,2,1)
for j=1:height(h_comb)
    semilogy([1:8], end_comb(j,:),'Color', [0 0 0.1 0.05]);
    hold on
    ylabel('Intensity')
    xlabel('Time (hr)')
    ylim([-5 18000])
    %set(gca, 'YScale', 'log')
end 

% Plotting normalized intensity for each droplet
subplot(1,2,2)
for j=1:height(h_comb)
    semilogy([1:8], norm_end_comb(j,:),'Color', [0 0 0.1 0.05]);
    hold on
    ylabel('Intensity')
    xlabel('Time (hr)')
    ylim([-5 350])
    %set(gca, 'YScale', 'log')
end 

% Plot normalized intensity against diameter
figure;
scatter(h_comb(h_comb >= 0)*3.24, end_comb(h_comb >= 0,8),'MarkerFaceColor', '#0073e6','MarkerFaceAlpha',0.2, 'MarkerEdgeAlpha',0);
hold on
dlm8 = fitlm(h_comb(h_comb >= 0)*3.24, end_comb(h_comb >= 0,8)) % Linear fit to show proportionality
plot([1:600], 107.67+[1:600]*28.148, 'r') %Plot linear fit
xlabel('Height of droplet (um)')
ylabel('Increase in OD')
title('Positive Controls')

% Plotting deviation from linear fit 
% Calculate percent error for each droplet
h_comb_new = h_comb(h_comb >= 0)*3.24;  % Filter by size
end_comb_new = end_comb(h_comb >= 0,8); % Filter by size
theo = h_comb_new*28.148 + 107.67; % Theoretical intensity from linear fit
perc_err = 100*abs(theo - end_comb_new)./theo; % Calculate percent error

% Calculate mean percent deviation over intervals
num_int = [0:50:600];
for i=1:length(num_int)
    mean_perc_err(i) = mean(perc_err(h_comb_new >= num_int(i) & h_comb_new <= num_int(i)+50));
end 

% Plot deviation form fit (mean and individual)
figure
scatter(h_comb_new, abs(perc_err),'MarkerFaceColor', '#0073e6','MarkerFaceAlpha',0.2, 'MarkerEdgeAlpha',0);
hold on
plot(num_int, mean_perc_err, 'Color','k')
%xline(100)
xlabel('Diameter (um)');
ylabel('Percent Error');
ylim([0 70])
