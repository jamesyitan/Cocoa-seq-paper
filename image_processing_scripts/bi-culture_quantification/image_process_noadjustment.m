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

    %to see individual photos, MAKE SURE TO CLICK ON NEW FIGURES AS
    %GENERATED OR CIRCLES WILL BE ON WRONG IMAGE
    %figure()
    %imshow(pc)
    %viscircles(centers_final, radii_final+5)
    
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

        [xgrid, ygrid] = meshgrid(1:size(mCherry,2), 1:size(mCherry,1));
        mask = ((xgrid-center(1,1)).^2 + (ygrid-center(1,2)).^2) <= radius.^2;
        intensity_droplet = mCherry(mask);

        intensity=sum(intensity_droplet);
        norm_intensity = intensity/(radii_final(n)^2*pi);
        %tracker_red_temp(n) = intensity;
        tracker_red_temp(n) = norm_intensity;
    end

    %add temporary trackers to the end of the previous interation
    tracker_green = cat(1,tracker_green,tracker_green_temp);
    tracker_red = cat(1,tracker_red,tracker_red_temp);

end

