function [cartvec bestJac polarvec index] = FindBestContrastVector(array1,array2,iter,res,image_tf,waitbar_tf)
%find best contrast vector in a couple of iterations 

%
%Example
%
%[cartvec bestJac polarvec index]=FindBestContrastVector(array1,array2,iter,res,imageyn)


%{


iter = 5;
res = 5
imageyn_tf = true;
waitbar_tf = true
%}

clear vector Jaccard Jacplot;



if waitbar_tf;
    hwb = waitbar(0,'iteration','Name','Find Best Contrast Vector');
end
    
measurements = res^2;

for k = 1:iter
    
    if waitbar_tf;
        l=1;
        waitbar((((k-1)*measurements)+l)/(measurements*iter),hwb,['iteration ',num2str(k),' of ',num2str(iter)]);
    end
    
    
    %create spherical unit vector in equal spaced positive directions
    if k == 1
        theta = linspace(0,pi,res);
        phi = linspace(0,pi,res);
    else
        rg = pi*0.5^k;  %range
        theta = linspace(direction(1,mr)-rg,direction(1,mr)+rg,res);
        phi = linspace(direction(2,mr)-rg,direction(2,mr)+rg,res);
    end
        
    direction = combvec(theta,phi);

    vector(1,:) = sin(direction(1,:)).*cos(direction(2,:));
    vector(2,:) = sin(direction(1,:)).*sin(direction(2,:));
    vector(3,:) = cos(direction(1,:));

    
    for l =1:measurements
        [Jaccard(l)] = sng_Rgb2BwContrast(array1,array2,vector(:,l),false);

        if waitbar_tf && floor(rem((((k-1)*measurements)+l),(measurements*iter/100))) == 0

            waitbar((((k-1)*measurements)+l)/(measurements*iter));
        end
        
        %(((k-1)*measurements)+l)/(measurements*iter)
    end

    [mx,mr]=max(Jaccard);
    
    index(k) = mr
    cartvec(1:3,k) = vector(:,mr)
    bestJac(k) = mx
    polarvec(1:2,k) = direction(1:2,mr)

    
    if image_tf
        Jacplot = reshape(Jaccard,size(theta,2),size(phi,2));
        figure;surf(theta,phi,Jacplot)
        xlabel('theta');ylabel('phi')
        xlim(sng_lim(theta))
        ylim(sng_lim(phi))
        colormap colorcube(200)
    end

end

if waitbar_tf;
    close(hwb);
end