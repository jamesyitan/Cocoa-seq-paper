%%Determine average value for "overrepresented" mCherry cells

mCherry = imread("microscopy_images/set_1/T=20_5-5_cy3.jpg");
mCherry = imbinarize(mCherry);
labeled_mCherry = bwlabel(mCherry);
blobProps_mCherry = regionprops(labeled_mCherry, 'Area','Centroid');
blobProps_mCherry_centroids = cat(1,blobProps_mCherry.Centroid);
blobProps_mCherry_area = cat(1,blobProps_mCherry.Area);

threshold_area = 50;

for i=1:length(blobProps_mCherry_area)
    if blobProps_mCherry_area(i) > threshold_area
           blobProps_mCherry_area(i) = NaN;
    end
end

blobProps_mCherry_radii = sqrt(blobProps_mCherry_area/pi);

imshow(mCherry)
hold on
viscircles(blobProps_mCherry_centroids,blobProps_mCherry_radii)
hold off

%%looks good? What area do we assign the overrepresented values with?

mean_blobProps_mCherry_area = mean(blobProps_mCherry_area(~isnan(blobProps_mCherry_area)));

desired_size = sqrt(13.4253);     % Adjust the desired area for replacement - using square to replace

%%
%%I suspect the fluorescence/cell is larger for mCherry cells than for GFP
%%cells (at least visually). Let's normalize that as well.

GFP = imread("microscopy_images/set_1/T=20_5-5_fitc.jpg");
GFP = imbinarize(GFP);
labeled_GFP = bwlabel(GFP);
blobProps_GFP = regionprops(labeled_GFP, 'Area','Centroid');
blobProps_GFP_centroids = cat(1,blobProps_GFP.Centroid);
blobProps_GFP_area = cat(1,blobProps_GFP.Area);

%using the same area threshold as above for mCherry, let's see if it gets rid of the large coalescing blobs 
for i=1:length(blobProps_GFP_area)
    if blobProps_GFP_area(i) > threshold_area
           blobProps_GFP_area(i) = NaN;
    end
end

blobProps_GFP_radii = sqrt(blobProps_GFP_area/pi);

imshow(GFP)
hold on
viscircles(blobProps_GFP_centroids,blobProps_GFP_radii)
hold off

%looks good, let's compare distributions (also need to run the section
%above too to get mCherry distribution)

%subsample
blobProps_GFP_area_subsample = datasample(blobProps_GFP_area(~isnan(blobProps_GFP_area)),400,'Replace',false);
blobProps_mCherry_area_subsample = datasample(blobProps_mCherry_area(~isnan(blobProps_mCherry_area)),400,'Replace',false);

%generate histogram to compare
figure()
histogram(blobProps_GFP_area_subsample,'FaceColor',[0.4660 0.6740 0.1880],'BinWidth',5)
xlabel('Area of single cells')
hold on
histogram(blobProps_mCherry_area_subsample,'FaceColor',[0.6350 0.0780 0.1840],'BinWidth',5)

avg_mCherry = mean(blobProps_mCherry_area(~isnan(blobProps_mCherry_area)));
avg_GFP = mean(blobProps_GFP_area(~isnan(blobProps_GFP_area)));

%%not very significant of a difference - no need to adjust

%% 

%%%test section to see what blob replacement looks like

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

imshow(mCherry)
figure()
imshow(edited_mCherry)
