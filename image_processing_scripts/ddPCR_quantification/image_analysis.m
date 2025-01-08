%collect data on NTC droplets
[blank1,blank_centers,blank_radii] = image_process('images/(#10-45)/blank.tif','images/(#10-45)/blank_f.tif',15,20,0.92);
blank_avg = nanmean(blank1);

%collect data on dilution
[dilution,dilution_centers,dilution_radii] = image_process('images/(#10-45)/dilution-5.tif','images/(#10-45)/dilution-5_f.tif',15,20,0.90);
num_total = length(dilution);

%set threshold from estimations from NTC droplets
threshold = 1400;
pos_tracker = dilution>threshold;
pos_centers = dilution_centers(pos_tracker,:);
pos_radii = dilution_radii(pos_tracker);
pos = dilution(pos_tracker);
num_pos = length(pos);
lambda = num_pos/num_total;

% imshow(imread('images/(#10-45)/dilution-5.tif'))
% viscircles(dilution_centers,dilution_radii)

% visualize
figure()
imshow(imread('images/(#10-45)/dilution-5_f.tif')*10)
viscircles(dilution_centers,dilution_radii,'Color','b')
viscircles(pos_centers,pos_radii)

%1 um is 1.08 pixel