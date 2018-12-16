function test
p=mfilename('fullpath')
rootDir=fileparts(p)

mexDir = [rootDir filesep 'mex' filesep computer('arch') filesep]
end