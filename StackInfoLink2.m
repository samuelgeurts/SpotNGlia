 v function [stackinfo] = StackInfoLink2(imageinfo)
%creates structure with stack information
%Version StackInfoLink2 20170616
%   changed to open CorNextStack instead of nextstack (from imageinfo)
imageinfo([imageinfo.CorNextStack] == 2) = [];


fk = find([imageinfo.CorNextStack]);
lk = [fk(2:end)-1,numel(imageinfo)];
stacksize = lk-fk+1;
for k = 1: numel(fk)
        
    firstim = imageinfo(fk(k)).name;
    lastim = imageinfo(lk(k)).name;
    %select all number right before '.tif'
    firstn = char(regexp(firstim,'\d+.tif','match'));  
    %remove the '.tif' part
    firstn = strrep(firstn,'.tif','');
    %select name without index
    name = strrep(firstim, [firstn,'.tif'],'');
    %if empty it is probably the first image and '000' is used

    %do the same for the last image
    lastn = char(regexp(lastim,'\d+.tif','match'));
    lastn = strrep(lastn,'.tif','');
    
    for k2 = 1:numel(lastn) - numel(firstn)
        firstn = [firstn,'0'];
    end    
    
    FileName = [name,'(',firstn,'-',lastn,')'];

    %for k2 = selec(1):selec(end)
    %    batchname{k2,1} = FileName;
    %end
    
    stackinfo(k).stackname = FileName;
    tmp = num2cell(stacksize);[stackinfo.stacksize] = tmp{:};
    clear imagenames
    for k2 = 1:stacksize(k)
        imagenames{k2} = imageinfo(fk(k)+k2-1).name;
    end
    stackinfo(k).imagenames = imagenames; 
    
end     
    
% %probeelsel met while
% k = 1;
% k2 = 1;
% while k <= numel(imageinfo)
%     %select all number right before '.tif'
%     firstn = char(regexp(imageinfo(k).name,'\d+.tif','match'));
%     %remove the '.tif' part
%     firstn = char(regexp(firstn,'\d+','match'));
%     %if empty it is probably the first image and '000' is used
%     if isempty(firstn)
%         firstn = '000';
%     end
%     
%     kf = k;
%     k = k + 1;
%     
%     while k <= numel(imageinfo) && imageinfo(k).nextstack == 0
%         k = k + 1;
%     end
%     
%     %retrieve last stackimage name
%     lastn = char(regexp(imageinfo(k-1).name,'\d+.tif','match'));
%     lastn = char(regexp(lastn,'\d+','match')); 
%     
%     stackinfo(k2).stackname = [firstn, '-' lastn];
%     stackinfo(k2).stacksize = k-kf;
% 
%     k2=k2+1
% end
%     

    
%old stack function
% for k = 1: fishnumber(end)
%     selec = find(fishnumber == k);    
%     firstim = imageinfo(selec(1)).name;
%     lastim = imageinfo(selec(end)).name;
%     %select all number right before '.tif'
%     firstn = char(regexp(firstim,'\d+.tif','match'));
%     %remove the '.tif' part
%     firstn = char(regexp(firstn,'\d+','match'));
%     %if empty it is probably the first image and '000' is used
%     if isempty(firstn)
%         firstn = '000';
%     end    
%     %do the same for the last image
%     lastn = char(regexp(lastim,'\d+.tif','match'));
%     lastn = char(regexp(lastn,'\d+','match'));
%     FileName = [firstn, '-' lastn];
% 
%     %for k2 = selec(1):selec(end)
%     %    batchname{k2,1} = FileName;
%     %end
%     
%     stackinfo(k).stackname = FileName;
% end
