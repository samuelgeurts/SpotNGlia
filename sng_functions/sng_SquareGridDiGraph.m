function [edges] = sng_SquareGridDiGraph(s)
%this function generates nodes and edges from a rectangular grid (matrix)
%direction is from left to right and connectivety is -45,0,45 degrees
%also the values from the matrix input in terms of nodes are calculated
%   s:  size of 2 dimensional input matrix



%!!!fix the starting point

%{
ex
matrix = rand(100,100);
[indcom val] = sng_SquareGridDiGraph(matrix)

%}

%size
%s = size(matrix)
%indices
[CC,RR] = meshgrid(1:s(2),1:s(1));

%RR = int32(RR);
%CC = int32(CC);



    %{
figure;scatter(CC(:),RR(:));axis ij
%}

%% first and last row
RRfr = RR(1,1:end-1);
RRlr = RR(end,1:end-1);
CCfr = CC(1,1:end-1);
CClr = CC(end,1:end-1);
%edges of first row
subfr = [RRfr(:),CCfr(:)];
efr1 = subfr + repmat([0 1],[length(subfr),1]); %one step right
efr2 = subfr + repmat([1 1],[length(subfr),1]); %one step down and right
%edges of last row
sublr = [RRlr(:),CClr(:)];
elr1 = sublr + repmat([0 1],[length(sublr),1]); %one step right
elr2 = sublr + repmat([-1 1],[length(sublr),1]); %one step up and right
%sub to ind first row
indfr = sub2ind(s,subfr(:,1),subfr(:,2));
indefr1 = sub2ind(s,efr1(:,1),efr1(:,2));
indefr2 = sub2ind(s,efr2(:,1),efr2(:,2));
%sub to in last row
indlr = sub2ind(s,sublr(:,1),sublr(:,2));
indelr1 = sub2ind(s,elr1(:,1),elr1(:,2));
indelr2 = sub2ind(s,elr2(:,1),elr2(:,2));
%edges of first en last row

%% indices except first and last row and last colum
RRmr = RR(2:end-1,1:end-1);
CCmr = CC(2:end-1,1:end-1);
%indices middle rows
submr = [RRmr(:),CCmr(:)];
%edges
e1 = submr + repmat([-1,1],[length(submr),1]);
e2 = submr + repmat([0,1],[length(submr),1]);
e3 = submr + repmat([1,1],[length(submr),1]);
%indices of edges
ind = sub2ind(s,submr(:,1),submr(:,2));
inde1 = sub2ind([s(1), s(2)],e1(:,1),e1(:,2));
inde2 = sub2ind([s(1), s(2)],e2(:,1),e2(:,2));
inde3 = sub2ind([s(1), s(2)],e3(:,1),e3(:,2));

%complete edges
edges = [[indfr,indefr1];[indfr,indefr2];[indlr,indelr1];[indlr,indelr2];[ind,inde1];[ind,inde2];[ind,inde3]];
%all values

%values = matrix(edges(:,2));


%{
%the adjacentcy matrix doesn have to be determined as digraph generates a
spart matrix
adj = zeros(24)
for j = 1:size(indc,1)
    adj(indc(j,1),indc(j,2)) = val(j)
end
%this does the same as accumarray
%}

%{
%in addition the adjacentcy matrix or the graph cal be computer from the
node and edges

adj = accumarray([indcom;indflc],val,[24 24]) %adjacentcy matrix
adj = accumarray(indcom,val,[prod(s),prod(s)],[],0,true); %sparse


G3 = digraph([indcom(:,1)],[indcom(:,2)],val);

%}






