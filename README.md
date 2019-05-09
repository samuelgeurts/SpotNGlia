# SpotNGlia
A Matlab toolbox described in the following paper:

**Reverse genetic screen reveals that Il34 facilitates yolk sac macrophage distribution and seeding of the brain**  
Laura E. Kuil, Nynke Oosterhof, Samuël N. Geurts, Herma C. van der Linde, Erik Meijering, Tjakko J. van Ham  
Published: MAR 2019   DISEASE MODELS & MECHANISMS   Volume: 12   Issue: 3    doi: 10.1242/dmm.037762

With this tool automated counting of microglia cells inside the midbrain of 3dpf zebrafish multistack brightfield images is performed. The algorithm can be distinguished in 5 main steps which have to be applied sequential. First, preprocessing which contains image corrections and enhancements on single fish images also the image slices of a single fish are merged such that in focus data is preserved. Second, image registration is applied, where a fish image is transformed such that it is aligned with a template fish with known orientation. Third, segmentation of the brain is applied to find the brain boundaries as we only have to look for microglia (red spots) in the brain area. Finally, a spotdetection method based on wavelets is applied to count microglia. A menu shows the results and enables manual correction. A csv file with microglia numbers is generated for further analysis.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

This package requires [Matlab](https://nl.mathworks.com/) R2015b or later. 

It requires the following toolboxes:
* [Signal Processing Toolbox](https://nl.mathworks.com/products/signal.html)
* [Image Processing Toolbox](https://www.mathworks.com/products/image.html)
* [Statistics and Machine Learning Toolbox](https://nl.mathworks.com/products/statistics.html)
* [Curve Fitting Toolbox](https://nl.mathworks.com/products/curvefitting.html)

It uses the following third party toolboxes and functions included in the [library folder](https://github.com/samuelgeurts/SpotNGlia/tree/backtTo1.5.3/library)
* [Image Alignment Toolbox](https://sites.google.com/site/imagealignment/)
* [normxcorr2_general](https://nl.mathworks.com/matlabcentral/fileexchange/29005-generalized-normalized-cross-correlation?focused=5228526&tab=function)

### Installing

Download the latest stable build of SpotNGlia to a local folder

[Releases](https://github.com/samuelgeurts/SpotNGlia/releases/tag/v1.6.0)

Run Matlab and add the folders and subfolders to path

## Demonstration
Run the following commands in order into your MATLAB command window

### initializes SpotNGlia object 
Shows decision menu for selection of the folder containing fish images and the folder to save data. Use the folder *testdata* for a demo.
```
obj = SpotNGlia;
```
### Associates image to fish
Shows a popup menu to choose a sorting method for the image folder. The images are sorted based on Date(capture time) or based on Name. Based on correlation, it is determined if adjoined images belonging to the same fish or not. As images of a single fish are normally imaged one after each other, sorting on *Date* is most used.
```
obj = obj.SliceCombination
```
A popup text appears with explaines the possibility to change the image-to-fish association or exclude images for the core algorithm:

Check slice combination in ImageInfo and StackInfo.
Apply corrections in "imageinfo.CorNextStack".
* Set value "1" for new fish
* Set value "2" for removing image from stack.
* Set value "0" if slice belongs to previous slice

### check image-to-fish association 
The variables 'obj.ImageInfo' and 'obj.StackInfo' are opened to check if the image to fish association is correct. 
At the column *CorNextStack* in the variable 'obj.ImageInfo' changes of the numbers can be made corresponding the legend above. In the variable 'obj.StackInfo' the results of the image to fish association is visible. If changes are made to obj.ImageInfo, obj.StackInfo will be overwritten when proceeding.
```
openvar('obj.ImageInfo')
openvar('obj.StackInfo')
```
### runs core algorithms
Runs all core algorithms in order after it is confirmed that obj.ImageInfo contains the correct image to fish association. A time bar shows up, displaying subsequently *Preprocession*, *Extended dept of field*, *Registration*, *Brain Segmentation* and *Spot Detection*.
```
obj = obj.CompleteProgram
```
### check results and make corrections
The fishes are shown and braincorrections and spotcorrections can be made
```
obj = obj.CheckFish
```
A screenshot of the user interface which can be used to check and correct brain edge and spot locations.

![User Interface](https://samuelgeurts.github.io/SpotNGlia/UserInterface.png)

A csv file is generated in the chosen save folder containing all uncorrected and corrected spot data.

## other usefull functions

### individual algorithm functions which has to be processed in order
```
obj = obj.PreProcession
obj = obj.ExtendedDeptOfField
obj = obj.Registration
obj = obj.BrainSegmentation
obj = obj.SpotDetection
````
### supporting properties 
These properties can be set after the SpotNGlia object is created.

Prevents to show screen that alerts for check of stackinfo and imageinfo.
```
obj.ImageInfoChecked_TF = true
```
Presetting obj.Sorting *Date* or *Name* prevents to show option menu.
```
obj.Sorting = 'Date'
```
Sets delimiter to ';' or ',' in csv file dependend on operating system.
```
obj.Delimiter = ';' 
```

## Versioning

We use [SemVer](http://semver.org/) for versioning.

## License

SpotNGLia is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
