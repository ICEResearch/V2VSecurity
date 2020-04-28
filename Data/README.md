# Data:
Data was gathered on center street in Provo, on the block between University
and Freedom (with the City Center Temple and NuSkin). The test points were 
as follows: 
* Alpha: V2I Scenario

* Foxtrot: This was "Bob", who was right behind "Alice"
* Golf: Eve was immediately behind Bob
* Hotel: Eve was about 1/2 city block behind Bob
* India: Eve was about 1 city block behind Bob
* Juliet: Eve was traversing in the direction opposite of Bob

LIST OF TEST POINTS: 

## GPSData: 
    Contains .kml files showing the path of the cars. These can be 
                plotted using Google Earth.
    Generated in the 'acquire' function
    Used in V2VGraphableData.m and V2IGraphableData.m

## GraphPwelchedData
    (I am honestly not sure what this stuff is for. I will try to keep an 
        eye out as I organize this repo, and will update this when I know)

## MetaData: 
    Contains .txt files with MetaData of timestamps, number of blocks,
                duration, etc... 
    Generated in the 'acquire' function

## PwelchedData:
    Contains .mat files with processed data gathered on the road. Note that 
        each test point is pwelched in its entirety, but there is also a 
        version that only uses the middle 32 carriers, which are later used
        for the coded portion (it needed to be a power of 2).

## CorrectlyAveragedData
    Contains .mat files that were, presumably, averaged somehow (correctly).
        Not sure what that means. 
    Generated in GenerateAveragedData.m (Likely with a length of 41) for 'averaged'
        Generated in GenerateGraphableData.m for 'jagged'
    Used in threshold_method.m (and through that, codedV2V.m) in order to 
        generate equivocation as a function of distance (See Fig. 9)
    
    CorrectlyAverageData from GenerateAveragedData.m
