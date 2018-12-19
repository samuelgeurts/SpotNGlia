function RegNew = RemoveDoubleSpots(Reg1, Reg2, limit)
%{
         Reg1 = Regions1
         Reg2 = Regions0{k3+1}
%}
rc1 = reshape([Reg1.Centroid], 2, numel(Reg1))';
rc2 = reshape([Reg2.Centroid], 2, numel(Reg2))';
src1 = size(rc1, 1);
src2 = size(rc2, 1);

distance = zeros(src1, src2);
for k70 = 1:src1
    distance(k70, :) = sqrt((rc2(:, 1) - rc1(k70, 1)).^2+(rc2(:, 2) - rc1(k70, 2)).^2);
end
%indexrc2 gives for every rc1 spot the closest rc2 point
[mindist, indexrc2] = min(distance, [], 2);

%sort the minimal distance in accending order
[sm, I] = sort(mindist);

RegNew = Reg1;

I(sm <= limit); %coupled indices of rc1
indexrc2(I(sm <= limit)); %coupled indices of rc2

%change spots in RegNew from Reg1 to Reg2 that has a higher peakvalue in Reg2
for k71 = 1:src1
    if sm(k71) <= limit
        p1 = Reg1(I(k71)).MPPeak;
        p2 = Reg2(indexrc2(I(k71))).MPPeak;
        if p2 >= p1
            RegNew(I(k71)) = Reg2(indexrc2(I(k71)));
        else
            %RegNew(I(k71)) = Reg1(I(k71));
        end
        %Reg1(I(k71)).Centroid
        %Reg2(indexrc2(I(k71))).Centroid
    end
end

%add spots that only occur in Reg2
rc2unique = 1:src2;
rc2unique(indexrc2(I(sm <= limit))) = [];
for k3 = 1:numel(rc2unique)
    RegNew(end+1) = Reg2(rc2unique(k3));
end
end