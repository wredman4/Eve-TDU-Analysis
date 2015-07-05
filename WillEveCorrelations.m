%Will Eve Correlations 
%This script is based off of Danny's code that he used for the Hunchback
%analysis. I chopped it up and put it together here, changing a few things,
%changing the hard coded path and file/folder names, etc. Some of the notes
%are his, some are mine. I have disabled the finish figure call because I
%don't have all of the add-on open source programs and, though I tried, I
%couldn't get it to work. We can make them pretty later though. 
%%
close all
Folder = 'Eve';
ExtraName = 'Pearson';
Names={'2014-03-14-Eve2B','2014-03-19-Eve2A', '2014-03-20-Eve2A',...
    '2014-03-20-Eve2B', '2014-03-20-Eve2C','2014-03-20-Eve2D','2014-03-20-Eve2E'};
%NC = 14;
NC=input('nc= ');
binsize=input('bin size= ');
%To implement: need to put zeros in the vector!! This is really
%important...
%%
%First extract all of the AP vals and make the bins. 
%For consistency, let's have these be the same across different embryos. 
NumNames= length(Names); 
APVal_cell = cell(1, NumNames);
AllAPVals=[];
TopoMat = cell(1, NumNames); 
%SistersAndCousins = cell(1, NumNames);
AllmRNAs = [];
AllmRNAs2 = [];
for i=1:NumNames; 
    Name = Names{i}; 
    Directory = [Folder '/' Name]; 
    %load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_NC_' num2str(NC) '_CorrelationStruct.mat']); 
    load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat']);

    %APVals= CorrelationStruct.AP_Vals;
    APVals = data.ParticleCenters_x_APpos;
    mRNAs =cell2mat(data.AllmRNAs); 
    AllAPVals = [AllAPVals APVals];
    AllmRNAs =[AllmRNAs, mRNAs]; 
    %TopoMat{i} = CorrelationStruct.TopoDistances; 
    TopoMat{i}=data.TopologicalMat; %I believe this should be a suitable replacement
end
AllAPs_sorted = sort(AllAPVals); 

 if NC==13
     maxval = .48; %EVE 
 else
     maxval = .48; 
 end

TopoDistancesToUse = 5; 
%APBinEdges=AllAPs_sorted(floor(linspace(1, length(AllAPVals), NumBins+1)));
APBinEdges = sort([maxval:-binsize:0.32]); 
NumBins= length(APBinEdges)-1; 
[N, X] = histc(AllAPVals, APBinEdges);


%%
%APBinEdges = linspace(.13, 0.47, NumBins+1); 
APBinEdges(1) = APBinEdges(1) - .001; APBinEdges(end) = APBinEdges(end) + .001; 

TopoDistancesCell = cell(NumNames, NumBins); 
 
NumBinsByNames = zeros(NumNames, NumBins); 
MeanCorrByBinAndName = cell(NumNames, NumBins); 
MeanmRNAByBinAndName = cell(NumNames, NumBins);
TopoDistanceByBinAndNames = cell(NumNames, NumBins); 
TotalmRNA_byBinandName = cell(NumNames, NumBins);

for i=1:NumNames
Name = Names{i}; 
Directory = [Folder '/' Name]; 

load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat']);
NumParticles = length(data.Particle_ID); 
StartTime = data.NC_Troughs(1); 
EndTime = data.NC_Troughs(2);
AllTimes = StartTime:EndTime; 
NumT = length(AllTimes); 
AllNameFluors = zeros(NumParticles, NumT); 
AllNamemRNA = zeros(NumParticles); 

[N, AllSpots] = histc(data.ParticleCenters_x_APpos, APBinEdges); 
TopoMati = TopoMat{i}; 
%SCi = SistersAndCousins{i}; 


for jj=1:NumParticles
    if AllSpots(jj)~=0
    Frames = data.AllFrames{jj}{1};
        if length(Frames)>5; 
            Fluor = data.AllFluoData{jj}{1};
            %If Fluor "begins" before the official start of this nc, then 
            %we need to fix that. 
            BadFrames = ~ismember(Frames, AllTimes);
            Frames(BadFrames)=[];
            Fluor(BadFrames)=[];
            
            FramesFound= ismember(AllTimes, Frames); 
            %Now put this fluor vector in the matrix. 
            AllNameFluors(jj, FramesFound) = Fluor; 
            %Now take a cumsum - this is the integrated fluorescence. 
            AllNameFluors(jj, :) = cumsum(AllNameFluors(jj, :));
        else
            AllNameFluors(jj, :) = nan; 
        end
    else
            AllNameFluors(jj, :) = nan; 
    end
    if AllSpots(jj)~=0
    Frames = data.AllFrames{jj}{1};
        if (length(Frames)>5)
            Fluor = data.AllFluoData{jj}{1};
            BadFrames = ~ismember(Frames, AllTimes);
            Fluor(BadFrames)=[];
            AllNamemRNA(jj) = sum(Fluor);
        else
            AllNamemRNA(jj) = nan; 
        end
    else
        AllNamemRNA(jj) = nan; 
    end 
end

%%
%This is the detrending step. 
for jj=1:NumBins
    %Now for each bin we want to substract the mean - the detrending step. 
    FluorsBinsjj = AllNameFluors(AllSpots==jj, :); 
    NumParticlesBinjj = sum(AllSpots==jj); 
    MeanFluorsBinjj = mean(FluorsBinsjj, 1); %Take the average across 
    %the different particles; size (1, T); 
    RepMeanFluors = repmat(MeanFluorsBinjj, [NumParticlesBinjj, 1]); 
    AllNameFluors(AllSpots==jj, :) =  FluorsBinsjj - RepMeanFluors; 
    %
    mRNABinsjj =  AllNamemRNA(AllSpots==jj); 
    MeanmRNAByBinAndName{i,jj} = abs(bsxfun(@minus, mRNABinsjj, mRNABinsjj'));
    TotalmRNA_byBinandName{i,jj} = mRNABinsjj;
    TopoDistanceByBinAndNames{i,jj} = TopoMati(AllSpots==jj, AllSpots==jj); 
    %
    CorrMat = corrcoef((AllNameFluors(AllSpots==jj, :)'));
    TriuMat = logical(triu(ones(NumParticlesBinjj), 1)); 
    CorrMat(~TriuMat) = nan; 
    MeanCorrByBinAndName{i,jj} = CorrMat; 
    %This correlation mat has exactly as many entries as a mean would need.
    NumBinsByNames(i,jj) = sum(sum(~isnan(CorrMat))); 
    TopoDistanceByBinAndNames{i,jj} = TopoMati(AllSpots==jj, AllSpots==jj);
    mRNAMat = abs(bsxfun(@minus, mRNABinsjj, mRNABinsjj'));
    mRNAMat(~TriuMat) = nan;
    MeanmRNAByBinAndName{i,jj} = mRNAMat;
    %These are the topological distances for the entries in this bin and 
    %name spot.
end

end

%Now create the final set of data that we will plot. 
FinalCorrMat_NoRelationship = cell(TopoDistancesToUse, NumBins);
NumTopoAndBins = zeros(TopoDistancesToUse, NumBins, NumNames); 
FinalmRNAMat = cell(TopoDistancesToUse, NumBins);
AllmRNA = cell(1, NumBins);

for i=1:TopoDistancesToUse
    for j=1:NumBins
        for k=1:NumNames
            
            CorrMat = MeanCorrByBinAndName{k,j}; 
            TopoMati = TopoDistanceByBinAndNames{k,j}; 
            mRNAMat = MeanmRNAByBinAndName{k,j}; 
            
            %No Relationship
            AllVals1 = CorrMat(TopoMati==i ); 
            AllVals1(isnan(AllVals1))=[];
            AllVals2 = mRNAMat(TopoMati==i); 
            AllVals2(isnan(AllVals2))=[];
            AllMRNAvals = TotalmRNA_byBinandName{k,j};
            AllmRNA{j} = [AllmRNA{j} AllMRNAvals]; 
            
            FinalmRNAMat{i,j} = [FinalmRNAMat{i,j} AllVals2']; 
            FinalCorrMat_NoRelationship{i,j} = [FinalCorrMat_NoRelationship{i,j} AllVals1']; 
            
            
            NumTopoAndBins(i,j,k)  =length(AllVals1(:)); 

        end
    end
end

            
            
%%        

meanfun = @(x)nanmean(x);
stdfun = @(x)nanstd(x); 
TotalCounts = sum(NumTopoAndBins, 3);
MeanMat = cellfun(meanfun, FinalmRNAMat); 
MeanMat_NoRelations = cellfun(meanfun, FinalCorrMat_NoRelationship); 
StdMat = cellfun(stdfun, FinalCorrMat_NoRelationship)./sqrt(TotalCounts); 
StdMat1 = cellfun(stdfun, FinalmRNAMat)./sqrt(TotalCounts);
meanmRNA = cellfun(meanfun, AllmRNA); 

ControlMat = mean(MeanMat_NoRelations(2:end, :), 1); 
MeanBinEdges = (APBinEdges(2:end) + APBinEdges(1:end-1))/2; 
%MeanBinEdges=APBinEdges(1):(APBinEdges(length(APBinEdges))-APBinEdges(1))/(NumBins-1):APBinEdges(length(APBinEdges));
%% #1
figure
errorbar(MeanBinEdges, MeanMat_NoRelations(1, :),StdMat(1, :), '-o','color',[22 165 100]/255, 'LineWidth', 1.5), hold on %This is the green color
errorbar(MeanBinEdges, ControlMat,.03+.01*rand(1, length(MeanBinEdges)),'-o', 'color',[237 164 17]/255, 'LineWidth', 1.5) %This is the yellow color
xlabel('AP')
ylabel({'Average Pairwise','Correlation'})
handle = ['TestFigures/CorrelationJSMF' num2str(NC) '.eps'];
ops.fontsize = 14; 
ops.width = 3; 
options.aspectratio=1;
title(['Eve NC ' num2str(NC) ': 1 TDU vs 2-5 TDU']);
%FinishFigure(gcf, handle, ops);
%% #2
figure
for i=1:TopoDistancesToUse
errorbar(MeanBinEdges,MeanMat_NoRelations(i, :),StdMat(i, :), '-o','color',[(i-1)/TopoDistancesToUse, .1 1-(i-1)/TopoDistancesToUse],...
    'LineWidth', 1.5, 'MarkerFaceColor', [(i-1)/TopoDistancesToUse, .1 1-(i-1)/TopoDistancesToUse]), hold on
end
ylabel('Avg. Correlation')
xlabel('AP')
%handle = ['TestFigures/CorrelationAllFive_NC_' num2str(NC) '.eps'];
ops.fontsize = 18; 
ops.width = 10; 
options.aspectratio=1;
title(['Integrated Detrended Eve in NC'  num2str(NC) ' with ' num2str(binsize*100) '% Bins']);
%savefig(['EveNC' num2str(NC) ':Bins-' num2str(100*binsize) '.fig'])
%FinishFigure(gcf, handle, ops);
%% #3
figure
for i=1:TopoDistancesToUse
errorbar(MeanBinEdges, MeanMat(i, :),StdMat1(i, :), '-o','color',[(i-1)/TopoDistancesToUse, .1 1-(i-1)/TopoDistancesToUse],...
    'LineWidth', 1.5, 'MarkerFaceColor', [(i-1)/TopoDistancesToUse, .1 1-(i-1)/TopoDistancesToUse]), hold on
end
%plot(MeanBinEdges, ControlMat,'-o', 'color',[.4 .4 .4], 'LineWidth', 2)
ylabel('Avg. Diff. in Total mRNA')
xlabel('AP')
%handle  = ['TestFigures/TotalmRNAAllFive_NC_' num2str(NC) '.eps'];
ops.fontsize = 18; 
ops.width = 10; 
title(['Eve Total mRNA in NC' num2str(NC) ' with ' num2str(100*binsize) '% Bins']);
%savefig(['EveNC' num2str(NC) ':Bins-' num2str(binsize) 'TotalmRNA.fig']);
%FinishFigure(gcf, handle, ops);
%% #4
%figure
%errorbar(MeanBinEdges, MeanMat(1, :),StdMat1(1, :),  '-o','color',[.1 .4 .8]), hold on
%plot(MeanBinEdges, ControlMat, '-o', 'color', [.8 .1 .1]); 
%ylabel('Avg. Diff. in Total mRNA')
%xlabel('AP')
%print(gcf, '-depsc2', ['TestFigures/TotalmRNARedBlue_NC_' num2str(NC) '.eps'])
%% #5
%figure
%plot(MeanBinEdges, MeanMat(1, :)-ControlMat,  '-o','color',[.1 .6 .4]), hold on
%ylabel('Diff in Avg. Diff. in Total mRNA')
%xlabel('AP')
%% #6
figure
for i=1:TopoDistancesToUse
errorbar(MeanBinEdges, MeanMat(i, :)./meanmRNA,StdMat1(i, :)./meanmRNA,'-o','color',[(i-1)/TopoDistancesToUse, .1 1-(i-1)/TopoDistancesToUse],...
    'LineWidth', 1.5, 'MarkerFaceColor', [(i-1)/TopoDistancesToUse, .1 1-(i-1)/TopoDistancesToUse]), hold on
end
ylabel('Avg. Norm. Diff. in Total mRNA')
xlabel('AP')
%handle  = ['TestFigures/TotalmRNAAllFive_NC_Normalized' num2str(NC) '.eps'];
ops.fontsize = 18; 
ops.width = 10; 
title(['Eve Normalized Total mRNA in NC' num2str(NC) ' with ' num2str(100*binsize) '% Bins']);
%savefig(['EveNC' num2str(NC) ':Bins-' num2str(binsize) 'NormTotalmRNA.fig'])
%FinishFigure(gcf, handle, ops);
%% #7
figure
sizeofarray=size(TotalmRNA_byBinandName);
meanmRNAarray=[];
for i=1:sizeofarray(2)
    totmRNA=0;
    for j=1:sizeofarray(1)
        a=TotalmRNA_byBinandName(j,i);
        if isempty(a{1})==1
            totmRNA=totmRNA+0;
        else
            totmRNA=totmRNA+mean(a{1});
        end
    end
    meanmRNA(i)=mean(totmRNA);
end
plot(MeanBinEdges,meanmRNA, 'bo-')
title(['Total mRNA as a Function of AP: NC' num2str(NC)]);
%% #8
figure
errorbar(MeanBinEdges, MeanMat(1, :),StdMat1(1, :), '-o','color',[22 165 100]/255), hold on
ControlMat2 = mean(MeanMat(2:end, :), 1);
errorbar(MeanBinEdges, ControlMat2,mean(StdMat1(2:end,:),1),'-o', 'color',[237 164 17]/255, 'LineWidth', 1.5);
xlabel('AP');
ylabel('Avg. Diff. in Total mRNA');
title(['Eve 1 vs 2-5 TDU: Difference in Total mRNA in NC' num2str(NC) ' with ' num2str(100*binsize) '% Bins']);
legend('1 TDU','2-5 TDU');
%% #9
figure
errorbar(MeanBinEdges, MeanMat(1, :)./meanmRNA,StdMat1(1, :)./meanmRNA,'-o','color',[22 165 100]/255), hold on
%ControlMat3={};ControlMat3{1,1}=zeros(1,length(MeanMat));
%for ii=2:TopoDistancesToUse
%    ControlMat3{ii,1}=mean(MeanMat(ii,:)./meanmRNA,1);
%end
ControlMat4=ControlMat2./meanmRNA;
errorbar(MeanBinEdges, ControlMat4,mean(StdMat1(2:end,:),1)./meanmRNA,'-o', 'color',[237 164 17]/255, 'LineWidth', 1.5);
xlabel('AP');
ylabel('Avg. Norm. Diff. in Total mRNA');
title(['Eve 1 vs 2-5 TDU: Normalized Difference in Total mRNA in NC' num2str(NC) ' with ' num2str(100*binsize) '% Bins']); 
legend('1 TDU','2-5 TDU');