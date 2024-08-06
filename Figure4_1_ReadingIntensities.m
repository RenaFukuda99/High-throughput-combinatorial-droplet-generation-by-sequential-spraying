%% Read in droplet intensities
% Read in droplet intensities and export to file.
% Recommended: generate a new file for each droplet "type" (i.e. different 
% file for positive control, fluorescent calibration droplets, etc.) 
% Data should then be passed to relevant MATLAB script (ex. flle for 
% fluorescent calibration droplets should be used for 
% Figure4_2_FluorescentCalibration"). 
clc;clear;

% Get directories and get relevant file names
D=dir();
Dark = D(contains({D.name},'4X_Dark.tif'));
Blue_300 = D(contains({D.name},'Blue_300_4X.tif'));
Blue_30 = D(contains({D.name},'Blue_30_4X.tif'));
Red_30 = D(contains({D.name},'Red_30_4X.tif'));
Red_300 = D(contains({D.name},'Red_300_4X.tif'));

% Set up darkfield thresholds
thresh_dia = 16000; %manually set intensity threshold

% Set up save file
filename = 'my_filename.xlsx'; %Recommend: name after image folder (ex. Fluorescent Calibration)
positions = []; %vector of XY positions of interest

%% Obtain intensity information from images
% Loop through each XY position

for j=1:height(Dark)
    
    % Get relevant images
    B30_now = Blue_30(contains({Blue_30.name}, strcat('XY',num2str(positions(j),'%03.f'))));
    B300_now = Blue_300(contains({Blue_300.name}, strcat('XY',num2str(positions(j),'%03.f'))));
    R30_now = Red_30(contains({Red_30.name}, strcat('XY',num2str(positions(j),'%03.f'))));
    R300_now = Red_300(contains({Red_300.name}, strcat('XY',num2str(positions(j),'%03.f'))));
    Dark_now =  Dark(contains({Dark.name}, strcat('XY',num2str(positions(j),'%03.f'))));

    % Create preliminary mask through intensity thresholding
    im_dark = imread(strcat(Dark_now(1).folder, '\', Dark_now(1).name)); 
    im_dark = imresize(im_dark, 2);
    mask = imfill(im_dark >= thresh_dia, 'holes');
    mask = bwareaopen(mask, 10);
    mask_raw = mask;

   % Remove edges to get rid of partial droplets
    mask(1:50,:) = 0; 
    mask(:,1:50) = 0;
    mask(2000:2044,:) = 0;
    mask(:,2000:2048) = 0;

    % Filter preliminary mask by circularity and size
    droplets = regionprops(mask, 'Centroid', 'Area', 'Circularity','MinorAxisLength', 'MajorAxisLength');
    droplets = droplets([droplets.Circularity] > 0.99);
    droplets = droplets([droplets.Area] <= 100000);
    droplets = droplets([droplets.Area] >= 20);

    % Catalog circle radii and centers
    diameters_1=([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2;
    centers = reshape([droplets.Centroid], 2, width([droplets.Centroid])/2)';
    radii = diameters_1/2;

    % Visualize circles
    imshow(im_dark);
    hold on
    viscircles(centers,radii);

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

    % Get droplet information and labels
    labeled_mask = bwlabel(mask, 8);
    droplets = regionprops(labeled_mask, 'Area', 'Circularity','MinorAxisLength', 'MajorAxisLength','Centroid');
    diameters=([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2;

    % Make smaller circular mask
    [X Y] = size(im_dark);
    [columnsInImage rowsInImage] = meshgrid(1:Y, 1:X);
    radii_sm = 0.7*((([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2)/2);
    filler_mask = zeros(X,Y);

    for k=1:length(droplets)
        centroid_x(k) = droplets(k).Centroid(2);
        centroid_y(k) = droplets(k).Centroid(1);
        circlePixels = (rowsInImage - centroid_x(k)).^2 + (columnsInImage - centroid_y(k)).^2 <= (radii_sm(k)).^2;
        filler_mask(circlePixels ~= 0) = true;
    end
    
    % Label smaller mask
    lab_mask_small = filler_mask.*labeled_mask;

    % Read blue channel images
    im_B30 = imread(strcat(B30_now(1).folder, '\', B30_now(1).name));
    bkgd_B30 = bsxfun(@times, im_B30, cast(~mask_raw, 'like', im_B30));
    bkgd_B30_num=median(nonzeros(bkgd_B30), 'all');
    sub_B30 = im_B30-bkgd_B30_num;
    drop_B30 = regionprops(lab_mask_small, sub_B30,'MeanIntensity');
    BlueInt_30 = [drop_B30.MeanIntensity]';

    im_B300 = imread(strcat(B300_now(1).folder, '\', B300_now(1).name));
    bkgd_B300 = bsxfun(@times, im_B300, cast(~mask_raw, 'like', im_B300));
    bkgd_B300_num=median(nonzeros(bkgd_B300), 'all');
    sub_B300 = im_B300-bkgd_B300_num;
    drop_B300 = regionprops(lab_mask_small, sub_B300,'MeanIntensity');
    BlueInt_300 = [drop_B300.MeanIntensity]';

    % Read red channel images
    im_R30 = imread(strcat(R30_now(1).folder, '\', R30_now(1).name));
    bkgd_R30 = bsxfun(@times, im_R30, cast(~mask_raw, 'like', im_R30));
    bkgd_R30_num=median(nonzeros(bkgd_R30), 'all');
    sub_R30 = im_R30-bkgd_R30_num;
    drop_R30 = regionprops(lab_mask_small, sub_R30,'MeanIntensity');
    RedInt_30 = [drop_R30.MeanIntensity]';

    im_R300 = imread(strcat(R300_now(1).folder, '\', R300_now(1).name));
    bkgd_R300 = bsxfun(@times, im_R300, cast(~mask_raw, 'like', im_R300));
    bkgd_R300_num=median(nonzeros(bkgd_R300), 'all');
    sub_R300 = im_R300-bkgd_R300_num;
    drop_R300 = regionprops(lab_mask_small, sub_R300,'MeanIntensity');
    RedInt_300 = [drop_R300.MeanIntensity]';

    % Get darkfield intensity information
    dark_int=[];
    
    for i=1:12
        
        % Read darkfield image
        im = imread(strcat(Dark_now(i).folder, '\', Dark_now(i).name));
        im = imresize(im, 2);

        % Get droplets
        im(im >= 40000) = NaN; % Remove edges 
        im = imtranslate(im, [(-0.5*i) -0.2*i]);
        drop_dark = regionprops(lab_mask_small, im, 'MeanIntensity');
        
        % Get intensity
        dark_int(:,i)= [drop_dark.MeanIntensity]';
        
    end 

    T=table(diameters, BlueInt_30, BlueInt_300, RedInt_30, RedInt_300, dark_int);
    writetable(T,filename,'Sheet',(j));

    clearvars -except Dark Dia Blue_30 Blue_300 Red_30 Red_300 thresh_dia filename positions
    close all
end