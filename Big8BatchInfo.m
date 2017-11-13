

%%
%{
obj = obj.NewPath(1)

filename = obj.SavePath
infoname = obj.InfoName
sngname = obj.SaveName
%}

P{1} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Csf1a B1'
I{1} = 'INFO_Csf1a B1'
S{1} = 'SNG_Csf1a B1'

P{2} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Csf1b B1'
I{2} = 'INFO_Csf1b B1'
S{2} = 'SNG_Csf1b B1'

P{3} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Il 34 B1'
I{3} = 'INFO_Il 34 B1'
S{3} = 'SNG_Il 34 B1'

P{4} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/WT B1'
I{4} = 'INFO_WT B1'
S{4} = 'SNG_WT B1'

P{5} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-12-01-Slco2b1-3dpf'
I{5} = 'INFO_2016-12-01-Slco2b1-3dpf'
S{5} = 'SNG_2016-12-01-Slco2b1-3dpf'

P{6} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2017-01-19-f13a1a.1-loc56'
I{6} = 'INFO_2017-01-19-f13a1a'
S{6} = 'SNG_2017-01-19-f13a1a'

P{7} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2017-01-26-fabp11a-Havcr2-3dpf'
I{7} = 'INFO_2017-01-26-fabp11a-Havcr2-3dpf'
S{7} = 'SNG_2017-01-26-fabp11a-Havcr2-3dpf'

P{8} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2017-02-03-p2rx3a-hcst-3dpf'
I{8} = 'INFO_2017-02-03-p2rx3a-hcst-3dpf'
S{8} = 'SNG_2017-02-03-p2rx3a-hcst-3dpf'

P{9} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-21 ch25h/Ch25h B2'
I{9} = 'INFO_Ch25h B2'
S{9} = 'SNG_Ch25h B2'

P{10} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-21 ch25h/WT B2'
I{10} = 'INFO_WT B2'
S{10} = 'SNG_WT B2'

P{11} = '/Volumes/Seagate Expansion Drive/PublicationBatch/20171026_DKO_NR_3dpf'
I{11} = 'INFO_20171026_DKO_NR_3dpf'
S{11} = 'SNG_20171026_DKO_NR_3dpf'



k1=1

load([P{k1}, '/', S{k1}, '.mat'])

obj = obj.NewPath(1)

filename = obj.SavePath
infoname = obj.InfoName
sngname = obj.SaveName


AlignedFish = imread([P{k1}, '/', 'AlignedFish', '/', obj.StackInfo(k1).stackname, '.tif']);



sng_SaveCell2TiffStack(CorrectedFishTemp, [TempFolderName, '/', obj.StackInfo(k1).stackname, '.tif'])


mfilename('fullpath')
[folder, name, ext] = fileparts(which('obj'));

S = dbstack('-completenames');
S(1).file


%{
for k2 = 1:8
    
    load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
    %load([P{1},'/',I{1},'.mat'], 'SpotsDetected')
    
    nspots = []
    for k1 = 1:numel(SpotsDetected)
        nspots(k1, 1) = numel(SpotsDetected{k1});
    end
    
    Sheet = [{obj.StackInfo.stackname}', {obj.StackInfo.stacksize}', num2cell(nspots)]
    title = {obj.InfoName, 'images', 'nspots'}
    
    ds = cell2dataset([title; Sheet]);
    export(ds, 'file', [P{k2}, '/', I{k2}, '.csv'], 'delimiter', ',')
end
%}
%{
 %writes to text file
 fileID = fopen([P{1},'/',I{1},'.txt'],'w');
 %fprintf(fileID,'%s %u\n',obj.StackInfo.stackname, nspots);
 fprintf(fileID,'%s %s %s \n','filename','images','spots');
 for k = 1:12
 fprintf(fileID,'%s %u %u \n',obj.StackInfo(k).stackname, obj.StackInfo(k).stacksize, nspots(k));
 end
 fclose(fileID)
%}

%{
 does not work
 fid = fopen([P{1},'/',I{1},'.csv'], 'w') ;
  fprintf(fid, '%s,', Sh{1,1:end-1}) ;
  fprintf(fid, '%s\n', Sh{1,end}) ;
  fclose(fid) ;
  dlmwrite('test.csv', Sh(2:end,:), '-append') ;
%}

