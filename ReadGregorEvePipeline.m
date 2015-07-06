% Read Gregor Eve Data Pipeline
%%
DataSource=[];
DataSource.Folder='Eve';
%DataSource.Folder='HcbData';
%'2014-03-14-Eve2B','2014-03-19-Eve2A',
Names={ '2014-03-20-Eve2A',...
    '2014-03-20-Eve2B', '2014-03-20-Eve2C','2014-03-20-Eve2D','2014-03-20-Eve2E'};
%Names={'6-26v2','11-15','11-12','11-10'};%'7-13',
for i = 1:length(Names)
    DataSource.Name=Names{i};
    data=readgregorpipeline_DKW_HcbOneSpot_Main_v2(DataSource);
end
