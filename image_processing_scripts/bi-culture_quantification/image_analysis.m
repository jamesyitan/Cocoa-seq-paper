%images from microscope were processed in ImageJ
%For FITC images, color channels were split and only green kept to remove
%blue background from microscope and "Brightness and Contrast" settings set
%between 25 and 77
%For CY3 images, color channels were split and only red kept 
%"Brightness and Contrast" settings set between 25 and 77
%For brightfield images, we need to enhance the black border around the droplets
%by "Process" > "Binary" > "Convert to mask" > "Binary" and "Dilate" the
%black around the droplets to allow for better contrast

%run image_process function for each condition imaged
[data_gfp_5_5,data_mCherry_5_5] = image_process('microscopy_images',40,80,0.97);
%[data_gfp_5_5_noadj,data_mCherry_5_5_noadj] = image_process_noadjustment('microscopy_images',40,65,0.96);

%curate data and remove NaN values
data_gfp_5_5 = data_gfp_5_5(~isnan(data_gfp_5_5));
data_mCherry_5_5 = data_mCherry_5_5(~isnan(data_mCherry_5_5));
%data_gfp_5_5_noadj = data_gfp_5_5_noadj(~isnan(data_gfp_5_5_noadj));
%data_mCherry_5_5_noadj = data_mCherry_5_5_noadj(~isnan(data_mCherry_5_5_noadj));

relabund_gfp = data_gfp_5_5./(data_gfp_5_5 + data_mCherry_5_5);
relabund_mCherry = data_mCherry_5_5./(data_gfp_5_5 + data_mCherry_5_5);

%relabund_gfp_noadj = data_gfp_5_5_noadj./(data_gfp_5_5_noadj + data_mCherry_5_5_noadj);
%relabund_mCherry_noadj = data_mCherry_5_5_noadj./(data_gfp_5_5_noadj + data_mCherry_5_5_noadj);

%% plot results


figure()
subplot(1, 2, 1)
histogram(randsample(relabund_gfp,500),'FaceColor',[0 0.4470 0.7410],'BinWidth',0.20)
ylim([0, 500])
hold on
xlabel('B.subtilis relative abundance - microscopy')
seq_data_bs_relabun = readtable("seq_Bs_relabund.csv",'Delimiter', '\t');
subplot(1, 2, 2)
histogram(randsample(table2array(seq_data_bs_relabun(:,2)),500),'FaceColor',[0.6350 0.0780 0.1840],'BinWidth',0.20)
ylim([0, 500])
xlabel('B.subtilis relative abundance - sequencing')