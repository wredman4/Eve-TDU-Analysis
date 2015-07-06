%Read Gregor data Pipeline
%This is used to mine the important data structures from Hernan's data
%structures. Load the lin.mat and CompiledParticles and you should be good
%to go
%WTR: All my comments got deleted the other day. Rather annoying. So I'll
%add a few, but I get the gist of the code for the most part (all my
%comments will start with (and end with if multiple lines) WTR). :WTR
function data=readgregorpipeline_DKW_HcbOneSpot_Main_v2(DataSource)

Folder = DataSource.Folder; 
Name = DataSource.Name; 
Directory = [Folder '/' Name]; 
   

%This path should go to where you are storing data locally on your
%computer. 
Files =  dir(fullfile(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory], '*.mat')); %WTR: Changed this to fit my own path
for i=1:length(Files)
    load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '/' Files(i).name])
    %load([Files(i).name]);
end
%%
%%
%First we feed schnitzcells, ElapsedTime,endnc14 and CompiledParticles to
%read_gregor_data. This gives us the full data structure (data) and
%actnucind and toberemoved. 
data=readgregorpipeline_DKW_HcbTwoSpot(schnitzcells,ElapsedTime,CompiledParticles, Directory);
%%


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data]=readgregorpipeline_DKW_HcbTwoSpot(schnitzcells,ElapsedTime,CompiledParticles, Directory)
%READ_GREGOR_DATA: This function is the bulk of the readgregorpipeline.
%Essentially it creates the data structure cell, which we use to construct
%the movie and tell things that happen in each frame. The following notes
%are Madhav's


NumberOfParticles  =length(CompiledParticles); 
maxFrame = length(ElapsedTime)+1;
%0. 
Schnitzcells_Visible = zeros(length(schnitzcells), maxFrame); 
for i=1:length(schnitzcells)
    Schnitzcells_Visible(i, schnitzcells(i).frames) = 1; 
end

%1. Make a complete array of compiled particles. 
 
AllCompiledParticles = zeros(NumberOfParticles, maxFrame);
NCOfParticles = zeros(1, NumberOfParticles); 
for i=1:NumberOfParticles
    AllCompiledParticles(i, CompiledParticles(i).Frame) = CompiledParticles(i).Fluo;
    CP_Frames = CompiledParticles(i).Frame; 
    ParticleVisibleButNucleusIsNot = ~ismember(CP_Frames, schnitzcells(CompiledParticles(i).Nucleus).frames); 
    
    ToZero = CP_Frames(ParticleVisibleButNucleusIsNot ); 
    
    %AllCompiledParticles_Full(i,CompiledParticles(i).Frame) = CompiledParticles(i).Fluo;
    AllCompiledParticles(i, ToZero ) =0;   
    NCOfParticles(i) = CompiledParticles(i).nc;

end
    


%%
%Calculate the nuclear cycle for each particle and for every frame in the
%data set. 
% 

%This can be done with out the clicking - each Compiled particle has which
%frame it is in!

%Try to load t

%This path should go to where you are storing compiled data locally on your
%computer. 
try load(['/Users/danielwells/Documents/PostDoc/Projects/CollectiveBehaviorInDrosophilaEmbryos/Matlab/CompiledData_Folder/' Directory  '_Peaks.mat']);
    load(['/Users/danielwells/Documents/PostDoc/Projects/CollectiveBehaviorInDrosophilaEmbryos/Matlab/CompiledData_Folder/' Directory  '_Troughs.mat']);
catch err
%WTR: I'm a bit confused by the above commands. I don't know of ever having
%done anything with any .mat files dealing with Peaks and Troughs. I'll
%have to ask Danny.Also, we only defines these files at the end of the code. :WTR
    close all
    plot(mean(AllCompiledParticles)), hold on
    NumberPeaks = input('How many peaks are there? ');


    fprintf('Click the peaks corresponding to the NC you want to analyze\n')
    fprintf('And then click the Troughs surrounding them (2*n) \n');
    fprintf('For the Troughs, just click (2*n) times, there is some lag\n');
    %WTR: This last text prompt was added as there was some confusion
    %without it 
    [NC_Peaks_Loc, NC_Peaks_Mag] = ginput(NumberPeaks); 
    [NC_Troughs_All, NC_DividingLines] = ginput(2*NumberPeaks); 


    if NC_Troughs_All(end) > size(AllCompiledParticles, 2)
        NC_Troughs_All(end) = size(AllCompiledParticles, 2); 
    end

    NC_Troughs_All = round(NC_Troughs_All); 
    NC_Peaks_Loc = round(NC_Peaks_Loc); 

    close all
end
PossibleNuclearCycles = [14 13 12 11 10 9]; %If we are going to have a nuclear cycle, it is going to be in this order. (WTR: This is due to the way in which Eruc and Hernan imaged). 
NumberOfNuclearCycles = length(NC_Peaks_Loc); 
NuclearCyclesPresent = sort(PossibleNuclearCycles(1:NumberOfNuclearCycles)); 



data={};
for i=1:NumberOfNuclearCycles
    NC = NuclearCyclesPresent(i); 
    ParticlesInNC = NCOfParticles == NC;
    StartSpot = NC_Troughs_All(2*i-1); 
    EndSpot = NC_Troughs_All(2*i);
    PeakSpot = NC_Peaks_Loc(i); 
    
    FramesOfNC = StartSpot:EndSpot; 
    data.ncPresent = 1; 
    data.NC_Frames = FramesOfNC; 
    
    %Particles_Visible_nc = find([CompiledParticles.nc]==NuclearCyclesPresent(i));     
    Particles_Visible_nc = intersect(find(sum(AllCompiledParticles(:, FramesOfNC)>0, 2)>5), find(ParticlesInNC)); %If any particle is visible for fewer than five frames, we are not even going to include it!
    
    Nuclei_Visible_nc = [CompiledParticles(Particles_Visible_nc).Nucleus]; 

    [Nuclei_Visible_nc, UniqueParticleLocations, AllNuclei_nc] = unique(Nuclei_Visible_nc); 
    NumUniqueNuclei= length(Nuclei_Visible_nc); 
    
    AllParticles = Particles_Visible_nc; 
    Particles_Visible_nc = Particles_Visible_nc(UniqueParticleLocations);
    ParticlePairs = cell(1, NumUniqueNuclei); 
    %WTR: Particle Pairs allows us to determine lineage, namely sisters,
    %mothers, and cousins. :WTR 
    for ii=1:NumUniqueNuclei
        ParticlePairs{ii} = (AllParticles(AllNuclei_nc == ii))' ;
    end



    

    %WTR: What follows is really a better version of makemeasurements.m
    %from last summer. :WTR
    %Find all of the other nuclei.
    Schnitzcells_Visible_nc = find(sum(Schnitzcells_Visible(:, FramesOfNC)>0, 2)>5); 
    
    nuc = Nuclei_Visible_nc;
    cenx = zeros(1, length(nuc)); 
    ceny = zeros(1, length(nuc)); 
    
    %First, fill in the centers with nuclei associated with particles. 
    for ii=1:length(nuc)
%          [schnitzcells(nuc(ii)).cenx(schnitzcells(nuc(ii)).frames==i)]
%          [schnitzcells(nuc(ii)).ceny(schnitzcells(nuc(ii)).frames==i)]
        [val, Peakspot] = min(abs(schnitzcells(nuc(ii)).frames - PeakSpot)); 
        cenx(ii) = schnitzcells(nuc(ii)).cenx(Peakspot);
        ceny(ii) = schnitzcells(nuc(ii)).ceny(Peakspot); 
    end
    
    ParticleCenters_x = zeros(1, length(nuc)); 
    ParticleCenters_y = zeros(1, length(nuc));
    ParticleCenters_x_APPos = zeros(1, length(nuc));
    ParticleCenters_x_APpos = zeros(1, length(nuc));
    for ii = 1:length(nuc)
        ParticleSpot = Particles_Visible_nc(ii);
        [val, Peakspot] = min(abs(CompiledParticles(ParticleSpot).Frame - PeakSpot)); 
        ParticleCenters_x(ii) = CompiledParticles(ParticleSpot).xPos(Peakspot); 
        ParticleCenters_y(ii) = CompiledParticles(ParticleSpot).yPos(Peakspot); 
        ParticleCenters_x_APPos(ii) = CompiledParticles(ParticleSpot).APPos(Peakspot);
        ParticleCenters_x_APpos(ii) = CompiledParticles(ParticleSpot).APpos(Peakspot);
    end
    %For now, we are going to assume that twins have the same location - it
    %should only be off by a little bit, and saves a lot of work!
    

    %Then fill the in with everything else. There will be repeats here, but
    %we are going to deal with them below, using unique. 
    AllNucleiVisible = [nuc Schnitzcells_Visible_nc']; 
    
    SchnitzcellCenters_x_i = zeros(1, length(Schnitzcells_Visible_nc)); 
    SchnitzcellCenters_y_i = zeros(1, length(Schnitzcells_Visible_nc)); 
    
    for ii=1:length(Schnitzcells_Visible_nc)
        [val, Peakspot] = min(abs(schnitzcells(Schnitzcells_Visible_nc(ii)).frames - PeakSpot)); 
        SchnitzcellCenters_x_i(ii) = schnitzcells(Schnitzcells_Visible_nc(ii)).cenx(Peakspot); 
        SchnitzcellCenters_y_i(ii) = schnitzcells(Schnitzcells_Visible_nc(ii)).ceny(Peakspot); 
        
    end
    
    %The ordering is really important here, because we are going to call
    %unique next, on the rows!
    TotalCenters = [cenx SchnitzcellCenters_x_i; ceny SchnitzcellCenters_y_i]';
    [UniqueRows, UniquePlaces, ~] = unique(TotalCenters, 'rows'); 
    

    
    nuc = AllNucleiVisible(UniquePlaces); 
    cenx = UniqueRows(:, 1); 
    ceny = UniqueRows(:, 2); 
    
    
    [SortNuc, NucPlaces] = sort(nuc); 
    NC = NuclearCyclesPresent(i); 
    data.nuc  = SortNuc; 
    data.NucleiCenters_x = cenx(NucPlaces)';
    data.NucleiCenters_y = ceny(NucPlaces)';
    NucleiCenters_x = cenx(NucPlaces)';
    NucleiCenters_y = ceny(NucPlaces)';
    [mem, locations] = ismember(Nuclei_Visible_nc, SortNuc); 
    MatchRawNucleiToParticles = locations; %These are indices that match the nuclei to individual particles. 
    data.MatchRawNucleiToParticles = locations;
    %data(t).nuc(data(t).WhichNucleiForParticles(1)) =
    %data(t).Nuclei_Visible_i(1). 
    %data.NucleiWithAssociatedParticles_ID = Nuclei_Visible_nc; %Which nuclei have a particle?
    data.Particle_ID=Particles_Visible_nc'; %This matches directly with Compiled Particles. 
    data.ParticleAssociatedNuclei_Centers_x = NucleiCenters_x(MatchRawNucleiToParticles);
    data.ParticleAssociatedNuclei_Centers_y = NucleiCenters_y(MatchRawNucleiToParticles);
    
    data.ParticleCenters_x = ParticleCenters_x; 
    data.ParticleCenters_y = ParticleCenters_y; 
    data.ParticleCenters_x_APPos = ParticleCenters_x_APPos; 
    data.ParticleCenters_x_APpos = ParticleCenters_x_APpos; 
     
    data.NC_Peak = NC_Peaks_Loc(i); 
    data.NC_Troughs = [NC_Troughs_All(2*i-1) NC_Troughs_All(2*i)]; 
    data.SisterPairs = ParticlePairs; 
    %data.AllCompiledParticles =cell2mat(data.SisterPairs);

    SisterIndexer = [];
    for ii=1:length(data.SisterPairs)
        SisterIndexer = [SisterIndexer ii*ones(1, length(data.SisterPairs{ii}))]; 
    end
    data.SisterIndexer = SisterIndexer; 
    %data.RawNucleiToAllParticles  =data.MatchRawNucleiToParticles(SisterIndexer); 
    %data.Nuclei_xLoc_OfAllParticles = data.ParticleAssociatedNuclei_Centers_x(SisterIndexer); 
    %data.Nuclei_yLoc_OfAllParticles = data.ParticleAssociatedNuclei_Centers_y(SisterIndexer); 
    
    TotalUniqueNucleiWithParticles = length(data.SisterPairs); 
    AllFluoData = cell(1, TotalUniqueNucleiWithParticles); 
    AllFrameData = cell(1, TotalUniqueNucleiWithParticles); 
    AllFinalmRNAs = cell(1,TotalUniqueNucleiWithParticles); 
    for ii=1:TotalUniqueNucleiWithParticles
        NumSisters = length(data.SisterPairs{ii}); 
        FluoProfiles =cell(1, NumSisters); 
        Frames = cell(1, NumSisters); 
        mRNAs = zeros(1,NumSisters); 
        for jj=1:NumSisters
            FluoProfiles{jj} = CompiledParticles(data.SisterPairs{ii}(jj)).Fluo; 
            Frames{jj} = CompiledParticles(data.SisterPairs{ii}(jj)).Frame; 
            try
            mRNAs(jj) = CompiledParticles(data.SisterPairs{ii}(jj)).TotalmRNA; 
            catch err
                1; 
            end
        end
        AllFluoData{ii} = FluoProfiles; 
        AllFrameData{ii} = Frames; 
        AllFinalmRNAs{ii} = mRNAs; 
    end
    
    data.AllFluoData = AllFluoData; 
    data.AllFrames = AllFrameData; 
    data.AllmRNAs = AllFinalmRNAs; 
    
    %Find sister nuclei 
    AllSisterNuclei = zeros(TotalUniqueNucleiWithParticles); 
    for ii=1:TotalUniqueNucleiWithParticles
        for jj=1:TotalUniqueNucleiWithParticles
            if ii~=jj
                if schnitzcells(data.nuc(MatchRawNucleiToParticles(ii))).P == schnitzcells(data.nuc(MatchRawNucleiToParticles(jj))).P
                    AllSisterNuclei(ii,jj) = 1; 
                    AllSisterNuclei(jj,ii) = 1;  
                end
            end
        end   
    end
   data.AllSisterNuclei = AllSisterNuclei; 
    
    %Make Parent/Child Matrix
    AllCousins = zeros(TotalUniqueNucleiWithParticles); 
    for ii=1:TotalUniqueNucleiWithParticles
        for jj=1:TotalUniqueNucleiWithParticles
            if ii~=jj
                %Then jj is the parent of ii
                %WTR:schnitzcells(schnitzcells(data.nuc(MatchRawNucleiToParticles(ii))).P).P
                %WTR:schnitzcells(schnitzcells(data.nuc(MatchRawNucleiToParticles(jj))).P).P
                try
                if ~isempty(schnitzcells(schnitzcells(data.nuc(MatchRawNucleiToParticles(ii))).P)) && ~isempty(schnitzcells(schnitzcells(data.nuc(MatchRawNucleiToParticles(jj))).P))
                
                    if schnitzcells(schnitzcells(data.nuc(MatchRawNucleiToParticles(ii))).P).P == schnitzcells(schnitzcells(data.nuc(MatchRawNucleiToParticles(jj))).P).P
                        AllCousins(ii,jj) = 1;
                        AllCousins(jj,ii) = 1;

                    end
                end
                catch err
                    1;
                end
            end
        end   
    end
    data.AllCousins = AllCousins; 
    
%     count=0; 
%     p = AllParents; 
%     pcell ={}; 
%     while ~isempty(all(p(:)))
%         pcell(count+1) = p; 
%         p = p*AllParents; 
%         count= count; 
%     end
%         
%     
%     %Find hamming distance between nuclei 
%     HammingDistance = zeros(TotalUniqueNucleiWithParticles); 
%     for ii=1:TotalUniqueNucleiWithParticles
%         for jj=1:TotalUniqueNucleiWithParticles
%             if ii~=jj
%                 
% 
%             end
%         end   
%     end
%     data.AllSisterNuclei = AllSisterNuclei;
%     
    
    
    
    
    
    [Neighbors, box ] = VornoiCell(data);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save('test_Neighbors_2.mat', 'Neighbors'); 
    % save('test_data_2.mat', 'data')
    data.box = box; 
    TopoMat = TopologicalNeighborDistance(Neighbors); 
    data.TopologicalMat = TopoMat; 
    data.ParticleNucleiTopoDistances = TopoMat(MatchRawNucleiToParticles, MatchRawNucleiToParticles); 
    
    %Orbitals = Orbitals_DKW(data, Neighbors);
    %save('test_Orbitals_2.mat', 'Orbitals')
    %data.AllOrbitals = Orbitals;

    
    %This path should go to where you are storing compiled data locally on your
%computer. 
    save(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat'], 'data'); %WTR: Changes made to fit my own path
    save(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory  '\_Peaks.mat'], 'NC_Peaks_Loc');
    save(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory  '\_Troughs.mat'], 'NC_Troughs_All');

    1;
    close all
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VORONI CELL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Neighbors, box ] = VornoiCell( data ) %WTR: This is like the old code as well. Should be fairly time intensive.  
%VornoiCell Calculates the Vornoi cells 


storage_cell=[];
box = [(min(data.NucleiCenters_x)-1) (max(data.NucleiCenters_x)+1);(min(data.NucleiCenters_y)-1) (max(data.NucleiCenters_y)+1)]; 

%[vx,vy]=voronoi([(data(bb).NucleiCenters_x)'],[(data(bb).NucleiCenters_y)']);%we need to define...
[v,c]=voronoin([(data.NucleiCenters_x)' (data.NucleiCenters_y)']);%we need to define...

  
if(~isempty(v))
    BorderCells = zeros(1, length(c)); 
    for ii=1:length(c); 
        vii = v(c{ii}, :); 

        if any((vii(:, 1)<box(1, 1)) + (vii(:, 1)>box(1, 2)) + (vii(:, 2)<box(2, 1)) + (vii(:, 2)>box(2, 2)))
            BorderCells(ii) = 1; 
        end
    end

    for ii=1:length(c)
        vlist=c{ii};
        storage_cell(ii).nverts=[vlist];
        for zz=1:length(vlist)
            storage_cell(ii).vertposx(zz)=v(vlist(zz),1);
            storage_cell(ii).vertposy(zz)=v(vlist(zz),2);
        end
        kk=1;
        for jj=1:length(c)
            vlisttemp=c{jj};
            %(length(intersect(vlist,vlisttemp))==2)


            shared_members = find(ismember(vlist, vlisttemp));

            vii = v(vlist(shared_members), :); 
            if BorderCells(ii) && BorderCells(jj)
                found = 0; 
                for iii=1:length(shared_members)
                    if (vii(iii, 1)> box(1, 1)) && (vii(iii, 1)< box(1, 2)) && ...
                            (vii(iii, 2)> box(2, 1)) && (vii(iii, 2)< box(2, 2))
                        found = 1; 
                    end
                end

                if found && (ii~=jj)
                        %(length(intersect(vlist,vlisttemp))==2)
                        storage_cell(ii).ncell(kk)=jj;
                        kk=kk+1;
                end
            else
            if length(shared_members) == 2
                [max_sm_1, spot1] = max(abs(v(vlist(shared_members(1)),:)));
                [max_sm_2, spot2] = max(abs(v(vlist(shared_members(2)),:)));

                max_sm_1 = ((max_sm_1)*(sign(v(vlist(shared_members(1)),spot1)))); 
                max_sm_2 = ((max_sm_2)*(sign(v(vlist(shared_members(1)),spot2))));  

                Condition = (((max_sm_1 > box(spot1, 1)) && (max_sm_1 < box(spot1, 2))) && ((max_sm_2 > box(spot2, 1)) && (max_sm_2 <box(spot2, 2))));

                if Condition
                        %(length(intersect(vlist,vlisttemp))==2)
                        storage_cell(ii).ncell(kk)=jj;
                        kk=kk+1;

                end
            end
            end
        end
    end
end
   
Neighbors = storage_cell; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Neighbor topological distances
function Matrix = TopologicalNeighborDistance(Neighbors)

NumNuclei = length(Neighbors); 
Matrix = zeros(NumNuclei); 
for i=1:NumNuclei
    numFound = 1; %(We always assume that the distance from a nucleus to itself is zero. 
    NeighborsFound = i;
    Spot = 1; 
    while (numFound < NumNuclei) && (Spot<25)
        NewPotentialNeighbors = unique([Neighbors(NeighborsFound).ncell]);
        NewNeighbors = NewPotentialNeighbors(~ismember(NewPotentialNeighbors, NeighborsFound)); 
        Matrix(i, NewNeighbors) = Spot; 
        NeighborsFound = [NeighborsFound NewNeighbors]; 
        numFound = numFound + length(NewNeighbors); 
        Spot = Spot + 1; 
    end
    if Spot==25
        NotFound = ~ismember(1:NumNuclei, NeighborsFound); 
        Matrix(i, NotFound) = 25; 
        
    end
end
end
        
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULATE ORBITALS%%%%%%%%%%%%%%%%%%%%
% 
% function Orbitals = Orbitals_DKW(data, Neighbors)
% %To remind myself: cell contains the information about individual nuclei.
% %Unclear at which time point. 
% 
% 
% %%%%%%%%%%Define parameters
% NumOrbitals = 10; 
% 
% Orbitals = cell(1, NumOrbitals); 
% 
% NumActivatedNuclei = length(data.NucleiCenters_x); 
% NeighborStorageCell = cell(1, NumActivatedNuclei); 
% PreviouslyConsideredNodes = cell(1, NumActivatedNuclei); 
% Nb = Neighbors; 
% 
% 
% for i = 1:NumActivatedNuclei
%     %Find the neighboring nuclei of the activated nuclei currently under qustion.
%     NeighborsOfNuclei = Nb(i).ncell;  
%     %Remove any of the neighbors that are not part of the activated list. 
% %       NeighborsOfNuclei(~ismember(NeighborsOfNuclei, AllActivatedNuclei)) = [];
%     NeighborStorageCell{i} = NeighborsOfNuclei; 
%     PreviouslyConsideredNodes{i} = [i NeighborsOfNuclei]; 
% end
% 
% Orbitals{1}=NeighborStorageCell; 
% 
% %Now we want to essentially repeat the procedure until we have found up to
% %the tenth neighbors, without repeat, of the neighbors involved. 
% 
% for i=2:10
%     NeighborCell = Orbitals{i-1}; 
%     NthStorageCell = cell(1, NumActivatedNuclei); 
%     for j=1:NumActivatedNuclei
% 
%         NthNeighbors = [Nb(NeighborCell{j}).ncell]; 
%         AlreadyConsidered  = PreviouslyConsideredNodes{j}; 
%         NthNeighbors(ismember(NthNeighbors, AlreadyConsidered)) = []; 
%         PreviouslyConsideredNodes{j} = [AlreadyConsidered NthNeighbors]; 
%         NthStorageCell{j} = NthNeighbors; 
%     end
%     Orbitals{i} = NthStorageCell; 
% end
% 
% end
% % 



%WTR: Below looks like a dead end. Some code that never got put into the
%pipeline? :WTR

function out = CondenseMultioutputCell(varargin)
out = [];
for i=1:length(varargin)
    out = [out varargin{i}]; 
end
out = unique(out); 
end













