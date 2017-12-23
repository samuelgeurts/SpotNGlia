




ZFParametersTemp = obj.ZFParameters(strcmp({obj.ZFParameters.stage}, 'SpotDetection'));
ZFParametersTemp = sng_zfinput(ZFParametersTemp, 0, 'SpotDetection', 'RgbToGray', 'ColorToGrayVector', [0; 1; 0], ''); %color selection
ZFParametersTemp = sng_zfinput(ZFParametersTemp, 0, 'SpotDetection', 'Wavelet', 'ScaleBase', 0.333, ''); %color selection
ZFParametersTemp = sng_zfinput(ZFParametersTemp, 0, 'SpotDetection', 'MultiProduct', 'MPlevels', [6:11], ''); %color selection
ZFParametersTemp = sng_zfinput(ZFParametersTemp, 0, 'SpotDetection', 'MultiProduct', 'MPthreshold',1000, ''); %color selection
ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MinSpotSize',25,''); %color selection
ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MaxSpotSize',410,''); %color selection
ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MinProbability',0.025,''); %color selection


obj.ZFParameters = ZFParametersTemp
obj = obj.SpotDetection
obj.CheckBrain



