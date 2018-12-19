            objb = SpotNGliaB
            objb = objb.getSNG
            objb = objb.LoadObjects

            quicksave(objb)
                        
            objb.all(@saveit)
            
            objb = objb.all(@NewPath,2:5,12)
            objb = objb.allpar(@SpotDetection2,[])
            objb = objb.all(@LoadAnnotations,5)
            
            objb = objb.all(@Dummy,1:5,12);          
            objb = objb.all(@CheckFish,5);          
            objb = objb.all(@SpotVal,1:5)
          
            for k1 =1:5
                [DS{k1}, lab{k1}] = DefineSpotThresholds(objb.obj(k1))
            end
            
        
            DSc = [DS{1};DS{2};DS{3};DS{4};DS{5}];
            labc = [lab{1};lab{2};lab{3};lab{4};lab{5}];
            a = prdataset(DSc(:, 1:8), labc);
            
            
            
            %% change hue, shift red to the middle
            DSc2 = DSc          
            DSc2(:,1) = DSc2(:,1) + 0.5 - (DSc2(:,1) >= 0.5)

            
           
%{
            DSc = DS{1};
            labc = lab{1};
            b = prdataset(DSc(:, 1:8), labc);
            b = b * scalem(b, 'variance'); %scaling
            wp = parzenc(b);

            
            DSc = [DS{2};DS{3};DS{4};DS{5}];
            labc = [lab{2};lab{3};lab{4};lab{5}];  
            c = prdataset(DSc(:, 1:8), labc);
            err = b * wp * testc %error of testset c on trained classifier wb
            err = c * wp * testc %error of testset c on trained classifier wb
%}            
            
            a = prdataset(DSc2(:, 1:8), labc);
            a = a * scalem(a, 'variance'); %scaling 
            
            [B,C] = gendat(a,0.5);         
            tic
            wp = dtc(B)
            toc
            
            err1 = B * wp * testc %error of testset c on trained classifier wb
            err2 = C * wp * testc %error of testset c on trained classifier wb

            save(['/Users/samuelgeurts/Desktop/','rbsvc5050', '.mat'], 'wp','a','B','C','err1','err2')          
            
            objb = objb.all(@CheckFish,5)    
            objb.all(@saveit,[])
            
            
            %%
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
    
            
            
            