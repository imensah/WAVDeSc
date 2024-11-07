clear
clc
close all


%load data for Noisy 300
DataUMI=readtable("/home/isabel/Desktop/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/Noisy_UMI_3000.csv");

%Convert the table into numeric
DataUMI=table2array(DataUMI);

%Transpose numeric data
DataUMI2 =DataUMI';


%Estimate the decompositon level
declev = [3000 2000];
declev_wf = 'db6';
declev_d = wmaxlev(declev,declev_wf);

%wavelet decomposition and reconstruction on noisy data using wavdec2
WAVDESC = wdenoise(DataUMI2,declev_d,"Wavelet","db6","DenoisingMethod","Bayes","ThresholdRule","Hard","NoiseEstimate","LevelDependent");
%WAVDESCrd = round(WAVDESC);


%Transpose data again
WAVDESCTT = WAVDESC';

%convert negative expressions to zero
%replace negative expressions by zero
xx = find(WAVDESCTT < 0);
xd = WAVDESCTT;
xd (xx ) = 0;
xy=find(xd < 0);

%save("/home/elisa/Desktop/Final_Glory/Final_simulation/Data/denoised_wavdesc/xd.csv")
dlmwrite('/home/isabel/Desktop/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/denoised_wavdesc_db6D.csv',xd)

%Compare noisy and truth datasets and find RMSE
compare3 = sum(DataUMI1(:) - xd(:));

RMSE11 = sqrt(mean((DataUMI1-xd).^2)); %= 0.6834
R2=sum(RMSE11);


%%%plot a heatmap and stackplot
%stackedplot(DataUMI1);
%figure;
%heatmap(xd);

%Plotting
noisy = DataUMI1(:,1);
denoised = xd(:,1);

figure(1) 
subplot(2,1,1) 
plot(noisy); 
title('Cell 1 noisy') 
xlabel 'Genes'; 
 
subplot(2,1,2) 
plot(denoised); 
title('Cell 1 denoised') 
xlabel 'Genes'; 
 

%%load data for Noisy 300RMSE12
TrueUMI=readtable("/home/elisa/Desktop/Final_Glory/Final_simulation/Data/UMI/1Truth_300.csv");

%Convert the true table into numeric
TrueUMI1=table2array(TrueUMI);


%Compare noisy and truth datasets and find RMSE
compare4= sum(TrueUMI1(:) - DataUMI1(:));
compare5= sum(TrueUMI1(:) - xd(:));


RMSE12 = sqrt(mean((TrueUMI1-xd).^2)); %= 0.6834
R3=sum(RMSE12);


%Plotting
True = TrueUMI1(:,1);
denoised = xd(:,1);

figure(2) 
subplot(2,1,1) 
plot(True); 
title('Cell 1 True') 
xlabel 'Genes'; 
 
subplot(2,1,2) 
plot(denoised); 
title('Cell 1 denoised') 
xlabel 'Genes'; 
