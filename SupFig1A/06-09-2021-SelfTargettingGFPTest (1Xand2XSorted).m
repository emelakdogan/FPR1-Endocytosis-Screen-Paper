%% Example script using CytometryExperiment
%% Read in the data
path='D:\Emel_Data\06-09-2021-SelfTargettingGFPTest (2XSorted)';
f=CytometryExperiment(path);
%% Display the sample labels
f.printSampleLabels;
%% Determine parameters to use for gating
f.printParamLabels;
%% Find SSC parameters
f.getParamNums('SSC-A');
%% Make a gate
f=f.setGateFromHeatmap;  % when given the choice, you can name the gate R1
%% Make another gate
%You can specify an input to the function such as f=f.setGateFromHeatmap('R1;R2'). 
%This will show you the scatter plot with gates R1 and R2 already applied 
%so that you can draw your new gate based on the cleaned up data.
f=f.setGateFromHeatmap;  % when given the choice, you can name the gate R2
%% View a gate overlaid on a heatmap 
%This function now shows the fraction of cells that lie within the gate.
%To apply sequential gates, merge their names with a semicolon in between ('R1;R2' would mean apply gates R1 and R2)
%You can give multiple gates as input. e.g.:f.showGatingHeatmap(2,'CD13ratio','R1');
%it will first apply gate R1 (removing dead cells), and then overlay gate 'CD13ratio' on top of the cleaned up heatmap.
f.showGatingHeatmap(15,'R2','R1'); 
%% Look at a scatter plots. Before compensation you can check how much leak on channel into the other, for ex ch 8 to ch 9 (PE and mChy)
f.heatmapGatedValues(1, 8, 11, 'R1'); % Sample number 1, parameters 3 and 4, gate R1
%% Compute a compensation value crosstalk mCh and APC
% Look first at the histograms of the channels: have a sense of where each
% of them falls into the other
% Find positive cells from the control sample to use to compute the
% crosstalk values: (for ex, for sample only APC (Ch 11), you will check the
% histogram of APC, see where the other channel that leaks (for ex mCh)
% falls on APC
sampleNum=3;
minVal=1000; maxVal=10000;  % range of intensities to use for the reference channel
ind=f.data0{sampleNum}(:,3)>minVal & f.data0{sampleNum}(:,3)<maxVal & f.getGatedIndices(sampleNum,'R1'); 
density_scatter_heatmap(f.data0{sampleNum}(ind,3),f.data0{sampleNum}(ind,4),100:100:10000,100:100:10000);
%% Set the compensation value:
f.compensation(4,3)=median(f.data0{sampleNum}(ind,4)./f.data0{sampleNum}(ind,3));
%% Apply the compensation:
f=f.correctValuesUsingCompensationFactors;
%% Check the result:
sampleNum=3;
f.heatmapGatedValues(sampleNum, 7, 11, 'R1');
%% If desired or needed, reset the compensation values (before trying an alternate compensation strategy)
f=f.resetCompensationValuesToIdentity;

%%
set(groot, 'defaultLineLineWidth', 2.0)
%% Look at a histogram for CD81 in SC495 with AF488-antiCD81
f.histogramGatedValues([13 14 31],7,'R1'); % Samples 1-5, parameter 3, gate R1
title('')
xlabel('FITC (AF488-antiCD81)')
%% Look at a histogram for FPR1 endocytosis at high fMLF concentrations- SC495 cells
f.histogramGatedValues([15:19 25],8,'R1'); % Samples 1-5, parameter 3, gate R1
title('')
xlabel('PE')
%% Print a summary table
%f.printResultTableOnePopulation('R1',8,'mean',11,'mean',8,'frac-2000','CD13ratio','gated');
%That will give you mean of channel 8, mean of channel 11, 
%the fraction of cells with at least 2000 in channel 8, 
%and the fraction of cells falling into the CD13ratio gate. 
%(all after first applying the R1 and R2 gates)
f.printResultTableOnePopulation('R1',8,'mean',8,'frac-19810');  % Use gate R1 as the primary gate, display the mean for channel 3,
                                                                           % fraction greater than 200 for channel 3, and the fraction within
                                                                           % R1 that are also within R2
%% save processed data
gates=f.gates; % Also saving the gates separately as a backup
%compensation=f.compensation; %backup
save ([path filesep 'processed data 1-29-2021.mat'],'f','gates');
%% load processed data
load ([path filesep 'processed data 1-29-2021.mat'],'f','gates'); %,'compensation'
