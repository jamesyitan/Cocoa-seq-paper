function [tracker, centers_final, radii_final] = image_process(file_bf,file_f,low_dim,high_dim,sens)

%%%This function intakes a directory name (condition) and takes all of the
%%%replicates (defined as a directory with one phase contrast, GFP
%%%fluorescence, and mCherry fluorescence each) and records the signal for
%%%all droplets from all of them. The images have all been processed and
%%%standardized in ImageJ to adjust the brightness and convert all the
%%%files to .jpgs.

bf = imread(file_bf);
f = imread(file_f);

[centers_final,radii_final] = imfindcircles(bf,[low_dim high_dim],'Sensitivity',sens,'ObjectPolarity','bright');
num_circles_final = length(centers_final);

% remove imcomplete circles that are on the boundaries
radius = high_dim;
% make sure dimensions match the image files you are using
for i = 1:num_circles_final
if (centers_final(i,1) < (radius)) || (centers_final(i,1) > (1344-radius)) 
    centers_final(i,:) = [NaN,NaN];
    radii_final(i) = NaN;
end
end

for i = 1:num_circles_final
    if (centers_final(i,2) < (radius)) || (centers_final(i,2) > (1100-radius))
        centers_final(i,:) = [NaN,NaN];
        radii_final(i) = NaN;
    end
end

%     %to see individual photos, MAKE SURE TO CLICK ON NEW FIGURES AS
%     %GENERATED OR CIRCLES WILL BE ON WRONG IMAGE
%     figure
%     imshow(pc)
%     viscircles(centers_final, radii_final)

%set up trackers that get reset for each replicate directory
tracker = zeros(num_circles_final,1);

%go through each droplet and calculate and track brightness
for n=1:num_circles_final
center = centers_final(n,:);
radius = radii_final(n);

[xgrid, ygrid] = meshgrid(1:size(f,2), 1:size(f,1));
mask = ((xgrid-center(1,1)).^2 + (ygrid-center(1,2)).^2) <= radius.^2;
intensity_droplet = f(mask);

intensity=sum(intensity_droplet);
norm_intensity = intensity/(radii_final(n)^2*pi);
tracker(n) = norm_intensity;
end

