function [tracker_green, tracker_red] = image_process(folder_name,low_dim,high_dim,sens)

%%%This function intakes a directory name (condition) and takes all of the
%%%replicates (defined as a directory with one phase contrast, GFP
%%%fluorescence, and mCherry fluorescence each) and records the signal for
%%%all droplets from all of them. The images have all been processed and
%%%standardized in ImageJ to adjust the brightness and convert all the
%%%files to .jpgs.

%%read in the replicate sets for one condition through directory names - IMPORTANT: "set_*"
%%must be in the directory name to be recognized
rep_files = dir([folder_name,'/','set_*']);
rep_files_cell = struct2cell(rep_files);
rep_filesname = string(rep_files_cell(1,:));
num_reps = length(rep_filesname);

%%set up universal trackers for data
tracker_green = zeros(0,0);
tracker_red = zeros(0,0);

%%loop for each replicate by going through each directory
for m=1:num_reps
    set_files = dir([folder_name,'/',convertStringsToChars(rep_filesname(m)),'/','*.jpg']);
    set_files_cell = struct2cell(set_files);
    set_filesname = string(set_files_cell(1,:));

    %order is in alphabetical, so bf first, then cy3, then fitc
    pc = imread([folder_name,'/',convertStringsToChars(rep_filesname(m)),'/',convertStringsToChars(set_filesname(1))]);
    mCherry = imread([folder_name,'/',convertStringsToChars(rep_filesname(m)),'/',convertStringsToChars(set_filesname(2))]);
    gfp = imread([folder_name,'/',convertStringsToChars(rep_filesname(m)),'/',convertStringsToChars(set_filesname(3))]);
    
    [centers_final,radii_final] = imfindcircles(pc,[low_dim high_dim],'Sensitivity',sens,'ObjectPolarity','bright');
    num_circles_final = length(centers_final);

    % remove imcomplete circles that are on the boundaries
    radius = low_dim;
    % make sure dimensions match the image files you are using
    for i = 1:num_circles_final
        if (centers_final(i,1) < (radius)) || (centers_final(i,1) > (1360-radius)) 
            centers_final(i,:) = [NaN,NaN];
            radii_final(i) = NaN;
        end
    end

    for i = 1:num_circles_final
            if (centers_final(i,2) < (radius)) || (centers_final(i,2) > (1024-radius))
                centers_final(i,:) = [NaN,NaN];
                radii_final(i) = NaN;
            end
    end

    %convert grayscale to binary black and white
    mCherry = imbinarize(mCherry);
    gfp = imbinarize(gfp);

    %Go through the mCherry image to downsize blobs on the peripherals of
    %the droplets that are getting overrepresented

    labeled_mCherry = bwlabel(mCherry);
    blobProps_mCherry = regionprops(labeled_mCherry, 'Area','Centroid');

    threshold_area = 50;
    desired_size = sqrt(13.4253);     % Adjust the desired area for replacement to what you determine in "blob_detector_test.m" - using square to replace
    edited_mCherry = mCherry;  % Create a copy of the binary image
    
    for i = 1:numel(blobProps_mCherry)
        if blobProps_mCherry(i).Area > threshold_area
            % Create a binary mask for the region to remove the large blob
            [row, col] = find(labeled_mCherry == i);
            minRow = min(row);
            maxRow = max(row);
            minCol = min(col);
            maxCol = max(col);
            remove_mask = false(size(edited_mCherry));
            remove_mask(minRow:maxRow, minCol:maxCol) = true;
            edited_mCherry(remove_mask) = false;
            % Create a binary mask for the region to add the replacement blob
            midRow = (minRow+maxRow)/2;
            midCol = (minCol+maxCol)/2;
            add_mask = false(size(edited_mCherry));
            add_mask(floor(midRow-desired_size/2):floor(midRow+desired_size/2), floor(midCol-desired_size/2):floor(midCol+desired_size/2)) = true;
            edited_mCherry(add_mask) = true;
        end
    end
    
    %%if you want to see the changes
    %imshow(mCherry)
    %figure()
    %imshow(edited_mCherry)

    %to see individual photos, MAKE SURE TO CLICK ON NEW FIGURES AS
    %GENERATED OR CIRCLES WILL BE ON WRONG IMAGE
    figure()
    imshow(pc)
    viscircles(centers_final, radii_final+5)
    
    %set up trackers that get reset for each replicate directory
    tracker_green_temp = zeros(num_circles_final,1);
    tracker_red_temp = zeros(num_circles_final,1);

    %go through each droplet and calculate and track brightness
    %value for GFP
    %figure
    %imshow(gfp)
    %viscircles(centers_final, radii_final+5)

    for n=1:num_circles_final
        center = centers_final(n,:);
        radius = radii_final(n)+5;

        [xgrid, ygrid] = meshgrid(1:size(gfp,2), 1:size(gfp,1));
        mask = ((xgrid-center(1,1)).^2 + (ygrid-center(1,2)).^2) <= radius.^2;
        intensity_droplet = gfp(mask);

        intensity=sum(intensity_droplet);
        norm_intensity = intensity/(radii_final(n)^2*pi);
        %tracker_green_temp(n) = intensity;
        tracker_green_temp(n) = norm_intensity;
    end
    
    %figure
    %imshow(mCherry)
    %viscircles(centers_final, radii_final+5)
    
    %go through each droplet and calculate and track brightness
    %value for mCherry (USING EDITED)
    for n=1:num_circles_final
        center = centers_final(n,:);
        radius = radii_final(n)+5;

        [xgrid, ygrid] = meshgrid(1:size(edited_mCherry,2), 1:size(edited_mCherry,1));
        mask = ((xgrid-center(1,1)).^2 + (ygrid-center(1,2)).^2) <= radius.^2;
        intensity_droplet = edited_mCherry(mask);

        intensity=sum(intensity_droplet);
        norm_intensity = intensity/(radii_final(n)^2*pi);
        %tracker_red_temp(n) = intensity;
        tracker_red_temp(n) = norm_intensity;
    end

    %add temporary trackers to the end of the previous interation
    tracker_green = cat(1,tracker_green,tracker_green_temp);
    tracker_red = cat(1,tracker_red,tracker_red_temp);

end

