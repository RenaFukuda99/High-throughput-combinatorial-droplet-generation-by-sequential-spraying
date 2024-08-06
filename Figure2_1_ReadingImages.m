%% Reading in fluorescent calibration droplets (Figure 2)
% Reading in raw calibration droplet intensities, which were formed using
% solutions of known fluorophore concentration. 
clc;clear;

% Get directories and image names
D=dir();
im_list = D(contains({D.name},'3_4X')); 

% Manually set intensity threshold (change per image set)
thresh = 10000; 

% Set up save file
filename = 'AlexaFluor.xlsx'; %or "Fluorescein"
table_tot = table;

%% Reading images

for j=1:height(im_list)

    % Load in image and create preliminary mask using intensity threshold
    im_maskref = imread(strcat(im_list(j).name));
    mask = imfill(im_maskref >= thresh, 'holes');
    mask = bwareaopen(mask, 10);
    mask_raw = mask;
    % imshow(mask) 

    % Remove edges to remove partial droplets
    mask(1:30,:) = 0; 
    mask(:,1:30) = 0;
    mask(2000:2044,:) = 0;
    mask(:,2000:2048) = 0;

    % Filter mask by circularity and size
    droplets = regionprops(mask, 'Centroid', 'Area', 'Circularity','MinorAxisLength', 'MajorAxisLength');
    droplets = droplets([droplets.Circularity] > 0.99);
    droplets = droplets([droplets.Area] <= 100000);
    droplets = droplets([droplets.Area] >= 1);

    % Get droplet centers and radii
    diameters_1=([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2;
    centers = reshape([droplets.Centroid], 2, width([droplets.Centroid])/2)';
    radii = diameters_1/2;

    % Make circular mask to label droplets
    [X Y] = size(im_maskref);
    [columnsInImage rowsInImage] = meshgrid(1:Y, 1:X);
    mask_temp = zeros(size(im_maskref));

    for k=1:length(droplets)
        centroid_x(k) = droplets(k).Centroid(2);
        centroid_y(k) = droplets(k).Centroid(1);
        circlePixels = (rowsInImage - centroid_x(k)).^2 + (columnsInImage - centroid_y(k)).^2 <= radii(k).^2;
        mask_temp(circlePixels ~= 0) = true;
        mask = mask_temp;
    end

    % Get circle mask information
    labeled_mask = bwlabel(mask, 8);
    droplets = regionprops(labeled_mask, 'Area', 'Circularity','MinorAxisLength', 'MajorAxisLength','Centroid');
    diameters=([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2;

    % Make smaller circular mask for intensity reading
    [X Y] = size(im_maskref);
    [columnsInImage rowsInImage] = meshgrid(1:Y, 1:X);
    radii_sm = 0.9*((([droplets.MajorAxisLength]'+[droplets.MinorAxisLength]')/2)/2);
    filler_mask = zeros(X,Y);

    for k=1:length(droplets)
        centroid_x(k) = droplets(k).Centroid(2);
        centroid_y(k) = droplets(k).Centroid(1);
        circlePixels = (rowsInImage - centroid_x(k)).^2 + (columnsInImage - centroid_y(k)).^2 <= (radii_sm(k)).^2;
        filler_mask(circlePixels ~= 0) = true;
    end
    
    % Label smaller labeled mask
    lab_mask_small = filler_mask.*labeled_mask;
    
    % Get image background
    bkgd_3 = bsxfun(@times, im_maskref, cast(~mask_raw, 'like', im_maskref));
    bkgd_3_num=median(nonzeros(bkgd_3), 'all');

    % Get mean intensity per droplet
    sub_3 = im_maskref-bkgd_30_num;
    drop_3 = regionprops(lab_mask_small, sub_3,'MeanIntensity');
    Int_3 = [drop_3.MeanIntensity]';
    
    im_num = repelem(j, height(drop_30));

    % Write to table and remove unneeded variables
    table_now=table(im_num', diameters, Int_3);
    table_tot = [table_tot; table_now];

    clearvars -except im_list Int_3 thresh filename table_tot bkgd_3_num
    close all
end

% Write intensities to file
writetable(table_tot,filename,'Sheet',extractAfter(im_list(1).folder, "594\") ); % If fluorescein, replace "594\" with "Fluorescein\" 
