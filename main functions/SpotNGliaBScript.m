            objb = SpotNGliaB
            objb = objb.getSNG
            objb = objb.LoadObjects

            quicksave(objb)
                        
            objb = objb.all('saveit')
            
            objb = objb.all(@NewPath,2:5,12)
            objb = objb.allpar(@SpotDetection2,[])
            objb = objb.all(@LoadAnnotations,[])
            
            objb = objb.all(@Dummy,1:5,12);         
            objb = objb.all(@SpotVal,3)
            objb = objb.all(@SpotVal,4)
            objb = objb.all(@CheckFish,4)
