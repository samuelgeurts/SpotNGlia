function [minC,C,Q] = sng_threshold2(Eh,Bh)
%% this function calculates the optimal threshold for 2 distributions,
% minimizes the error
% only works for RGb with max value 255

%C is the threshold coordinate
%Q is good/total classifiers  
% Eh=Eh1
% Bh=Bh1
    for i=1:size(Eh,1)
    %look which distribution is left or right and compute the
    %falsepostive function

        for j = 1:size(Eh,2)
              F1(i,j) = sum(Bh(i,1:j))+sum(Eh(i,j:end));
              F2(i,j) = sum(Eh(i,1:j))+sum(Bh(i,j:end));
        end

        C(i,:)=[F1(i,:) F2(i,:)];

        %find the ratio good/false positives    
        Q(i)=1-(min(C(i,:))/sum([Eh(i,:) Bh(i,:)]));

        %find minimum C coordinate for threshold 
        %mean is used to find a minimum when more values are found
        %round is used to get an integer
        minC(i)=round(mean(find(C(i,:)==min(C(i,:)))));    

        if minC(i) >= 255
            minC(i)=minC(i)-255;
        end
    end
end

