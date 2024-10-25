% Generate concentration map (x,y in um and concentrations)
concentrationMap = generateConcentrationMap(csvread('xy_25um_c1c2c3.csv',0,0));
mapcentroid = concentrationMapCentroid(concentrationMap);
imshow(concentrationMap)