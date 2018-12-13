
function [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching(SpotCom,SpotAnn,limit)

%compares 2 lists of 2d coordinates
%matches coordinates based on distance up to a given limit in pixels

%Example   [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching...
%                    (SpotCom,obj.SpotInfo(k1).AnnotatedSpots,10);       
                
%{
SpotCom = SpotCom{k1}
SpotAnn = obj.SpotInfo(k1).AnnotatedSpots
limit = 10
%}


    lpc = size(SpotCom,1);
    lpa = size(SpotAnn,1);
    distance = zeros(lpc,lpa);

    for k70 = 1:lpc
        distance(k70,:) = sqrt((SpotAnn(:,1)-SpotCom(k70,1)).^2 + (SpotAnn(:,2)-SpotCom(k70,2)).^2);
    end

    %find for every computed spot the closest annotated point
    [mindist,indexAnn] = min(distance,[],2);
    %sort the minimal distance in accending order
    [~,I] = sort(mindist);

        %   figure;bar(1:length(SpotCom),mindist(I));xlim([1,length(SpotCom)])

    % closest distance matrix
    %linkmatrix(sub2ind([lpc,lpa],[1:lpc]',indexAnn(I))) = 1

    %% select Correct, FalsePostives and NegativePositives
    k80=0;
    k81=0;
    link = ([I,indexAnn(I),mindist(I)]);

    Correct = zeros(lpc,2);
    Correctlinks = zeros(lpc,1);
    FalsePos = zeros(lpc,2);
    FalseNeg = zeros(lpc,2);
    
    for k71 = 1:lpc
        %'1' if link is unique compared to previous links
        if  (k71 == 1) | ~(link(k71,2) == link(1:k71-1,2));
            link(k71,4) = 1;
        else
            link(k71,4) = 0;
        end

        %'1' if distance is smaller than limit, else '0'
        if  link(k71,3) <= limit
            link(k71,5) = 1;
        else
            link(k71,5) = 0;
        end

        %select all correct found computed spot coordinates
        if (link(k71,4) == 1) && (link(k71,5) == 1)
            k80 = k80 + 1;
            Correct(k80,:) = SpotCom(link(k71,1),:);
            Correctlinks(k80) = link(k71,2);
        end
        %select all incorrect found spot coordinates
        if (link(k71,4) == 0) || (link(k71,5) == 0)
            k81 = k81 + 1;
            FalsePos(k81,:) = SpotCom(link(k71,1),:);

        end
    end

    Correct(k80+1:end,:) = [];
    Correctlinks(k80+1:end) = [];
    FalsePos(k81+1:end,:) = [];
    
    k82=0;    
    for k72 = 1:lpa
        %select all spot coordinates that are not found
        %Correction made at 20171007
        %if exist('Correctlinks','var') & ~(k72 == Correctlinks)
        if exist('Correctlinks','var') && ~ismember(k72,Correctlinks)            
            k82 = k82 + 1;
            FalseNeg(k82,:) = SpotAnn(k72,:);
        end        
    end

    FalseNeg(k82+1:end,:) = [];

end