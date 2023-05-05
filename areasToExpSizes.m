function expSizes = areasToExpSizes(data)
 areas = table2array(data(:,{'Area__m_2_'}));
 % - remove N/A from data
 areas = areas(~isnan(areas));
 % - sort the data. areas.(str) is a matlab stucture, like a wardrobe with different
 % - shelves, one for each (str)
 areas = sort(areas);
 expSizes = areaToVol(areas);
return

function s = areaToVol(a)
 s = ((4/3).*pi.*(sqrt(a/pi)).^3);
return
