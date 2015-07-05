%% Mutual Information verison 3
%   To address the problem of not being able to tell what part of the
%   mutual information values were from the geometry of the configuration
%   of Eve (stripe) and what were from being topologically close, this
%   script bins the nuclei and then calculates the mutual information
%   among all two pairs up to 5 TDU away in each bin. The values within a
%   given bin for a specific TDU are then averaged and plotted, in a
%   similar fashion to Danny's figures in his write-up. This script
%   generated figure 4 in the new write-up. 
%%
tic
Folder = 'Eve'; 
Names={'2014-03-14-Eve2B','2014-03-19-Eve2A', '2014-03-20-Eve2A',...
    '2014-03-20-Eve2B', '2014-03-20-Eve2C','2014-03-20-Eve2D',...
    '2014-03-20-Eve2E'}; 
NC=14;maxTDU=5;
NumNames=length(Names);
BinSize=0.015; %or 0.010. should roughly be equal to a nuclear width (nw)
Bins=0.32:BinSize:0.48;%these were determined by Bothma et al.'s definition of stripe
%% MI Calculation
%   The mutual information of each nuclei of each bin of each embryo is
%   calculated and then saved to MIbyNamebyBin. 
MIbyNamebyBin=struct('Bin',{});
for ii=1:NumNames
    
     Name = Names{ii}; 
     Directory = [Folder '/' Name];
     load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat'])
     APpos=data.ParticleCenters_x_APpos;
     
     dataMI=struct('Bin',[],'MI',{});
     dataBin=struct('MI',{});
     %%
     X=zeros(length(APpos),length(data.NC_Frames)); 
     for jj=1:length(APpos)
         frames=data.AllFrames{jj}{1};
         frames=1+frames-min(data.NC_Frames);
         frames(frames<1)=[];
         X(jj,frames)=1;%X stores when in time a given nucleus is on
     end
     %% 
     % I have since found a more efficient way of doing this, but for the
     % time being, since it works, I'll keep it. 
     for jj=1:length(Bins)-1
         partInbin=find(APpos>Bins(jj) & APpos<Bins(jj+1));% all particles in bin
         if length(partInbin)~=1 && isempty(partInbin)~=1        
             for kk=1:length(partInbin)            
                 TDUMIarray=zeros(1,maxTDU);                 
                 for ll=1:maxTDU                 
                     neighbors=find(data.ParticleNucleiTopoDistances(partInbin(kk),:)==ll);                     
                     if isempty(neighbors)~=1                     
                         sameTDUarray=zeros(1,length(neighbors));                        
                         for mm=1:length(neighbors)                       
                             if isempty(intersect(neighbors(mm),partInbin))~=1 %making sure the neighbor is in the bin
                                sameTDUarray(mm)=MI_v3(partInbin(kk),neighbors(mm),X);
                            end                         
                         end                         
                         TDUMIarray(ll)=mean(sameTDUarray);%averaging the MI of all neighbors for a given particle for a given TDU                    
                     else                         
                         TDUMIarray(ll)=nan;                     
                     end
                 end                 
                 dataMI(partInbin(kk)).MI=TDUMIarray;
                 dataMI(partInbin(kk)).Bin=jj;
             end
         end
     end     
     %% 
     for jj=1:length(Bins)-1     
         ind=find([dataMI.Bin]==jj);         
         if isempty(ind)~=1         
             binMIvalues=[dataMI(ind).MI];
             binMIvaluesarray=zeros(1,maxTDU);             
             for kk=1:maxTDU             
                 TDUarray=zeros(1,length(binMIvalues)/maxTDU);
                 numberofdatasets=1:length(binMIvalues)/maxTDU;
                 TDUarray(numberofdatasets)=binMIvalues(kk+(numberofdatasets-1)*maxTDU);
                 TDUarray=TDUarray(~isnan(TDUarray));                
                 if isempty(TDUarray)~=1                 
                     binMIvaluesarray(kk)=mean(TDUarray);                 
                 else                     
                     binMIvaluesarray(kk)=nan;                
                 end
             end             
             dataBin(jj).MI=binMIvaluesarray;
         else
             dataBin(jj).MI=[nan,nan,nan,nan,nan];
         end
     end     
     MIbyNamebyBin(ii).Bin=dataBin;
end
%% Averaging MI by Bin
AverageMI=zeros(maxTDU,length(Bins)-1);%These will be what we use to plot the data
STDmatrix=zeros(maxTDU,length(Bins)-1);
for ii=1:length(Bins)-1
    for jj=1:maxTDU
        array2=zeros(1,NumNames);
        for kk=1:NumNames
            mi=MIbyNamebyBin(kk).Bin(ii).MI;
            if isempty(mi)==1 || isnan(mi(jj))==1
                array2(kk)=nan;
            else
                array2(kk)=mi(jj);
            end
        end
        AverageMI(jj,ii)=mean(array2(~isnan(array2)));
        STDmatrix(jj,ii)=std(array2(~isnan(array2)));
    end
end
%% Plotting
close all
%%  MI plot all 5 TDU
figure
MeanBins= (Bins(2:end) + Bins(1:end-1))/2;
for ii=1:maxTDU %same colors and parameters as in Danny's figures
    errorbar(MeanBins,AverageMI(ii,:),STDmatrix(ii,:),'-o','color',[(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU],...
    'LineWidth', 1.5, 'MarkerFaceColor', [(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU]), hold on
end
xlabel('AP');
ylabel('Averaged Mutual Information');
title(['Eve Mutual Information Analysis with ' num2str(BinSize*100) '% Bins']);
legend('1-TDU','2-TDU','3-TDU','4-TDU','5-TDU');
%%  MI plot 1 TDU vs. 2-5 TDU
figure 
errorbar(MeanBins, AverageMI(1, :),STDmatrix(1,:), '-o','color',[22 165 100]/255), hold on
two2fiveTDU=mean(AverageMI(2:end, :), 1);
errorbar(MeanBins,two2fiveTDU,mean(STDmatrix(2:end,:),1),'-o','color',[237 164 17]/255);
xlabel('AP');
ylabel('Averaged Mutual Information')
title(['Averaged 2-5 TDU MI vs. 1 TDU with ' num2str(BinSize*100) '% Bins' ]);
legend('1-TDU','2-5 TDU');
%% Normalized MI plot
figure 
for ii=1:maxTDU
    errorbar(MeanBins,AverageMI(ii,:)./mean(AverageMI),STDmatrix(ii,:)./mean(AverageMI),...
        '-o','color',[(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU],'LineWidth', 1.5,...
        'MarkerFaceColor', [(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU]), hold on
end
xlabel('AP');
ylabel('Normalized Averaged Mutual Information');
title(['Eve Normalized Mutual Information Analysis with ' num2str(BinSize*100) '% Bins']);
legend('1-TDU','2-TDU','3-TDU','4-TDU','5-TDU');
%% Mean Time On plot
figure
meantime=zeros(1,length(Bins)-1);
for jj=1:length(Bins)-1
    partInbin=find(APpos>Bins(jj) & APpos<Bins(jj+1)); 
    X2=X(partInbin,:);
    meantime(jj)=mean(mean(X2'));
end
plot(MeanBins,meantime,'ko-'); hold on


toc
            
         
                 