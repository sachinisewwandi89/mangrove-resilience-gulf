// 01_annual_composites — Phase 1: annual dry-season NBR series

// ---- 1. Study area ----
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -12.5]);

// ---- Mangrove extent guide layer ----
var gmw = ee.ImageCollection('projects/sat-io/open-datasets/GMW/annual-extent/GMW_MNG_2020')
            .mosaic()
            .clip(gulf)
            .selfMask();
Map.addLayer(gmw, {palette: ['00FF00'], opacity: 0.8}, 'Mangrove extent 2020');

// ---- 2. Cloud mask for Landsat Collection 2 (QA_PIXEL bits) ----
function maskL2(img) {
  var qa = img.select('QA_PIXEL');
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);      // cloud
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);     // cloud shadow
  var dilated = qa.bitwiseAnd(1 << 1).neq(0);    // dilated cloud
  var mask = cloud.or(shadow).or(dilated).not();
  return img.updateMask(mask);
}

// ---- 3. Scaling + NBR, per sensor ----
// NBR = (NIR - SWIR2) / (NIR + SWIR2)
// L5/L7: NIR=SR_B4, SWIR2=SR_B7 | L8/L9: NIR=SR_B5, SWIR2=SR_B7
function prepTMETM(img) {   // Landsat 5 & 7
  var sr = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  var nbr = sr.normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR');
  return nbr.copyProperties(img, ['system:time_start']);
}
function prepOLI(img) {     // Landsat 8 & 9
  var sr = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  var nbr = sr.normalizedDifference(['SR_B5', 'SR_B7']).rename('NBR');
  return nbr.copyProperties(img, ['system:time_start']);
}

// ---- 4. Merge all sensors into one NBR collection ----
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(gulf).map(maskL2).map(prepOLI);
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(gulf).map(maskL2).map(prepOLI);

var allNBR = l5.merge(l7).merge(l8).merge(l9);

// ---- 5. One dry-season composite per year (May–Oct) ----
var years = ee.List.sequence(1987, 2025);
var annual = ee.ImageCollection.fromImages(
  years.map(function(y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    return allNBR.filterDate(start, end)
      .median()
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);

// ---- 6. Sensitivity check A: view one composite (2016) ----
var img2016 = annual.filter(ee.Filter.eq('year', 2016)).first();
Map.setCenter(135.70, -15.10, 11);
Map.addLayer(ee.Image(img2016), 
  {min: -0.3, max: 0.7, palette: ['8c510a','f6e8c3','5ab4ac','01665e']},
  'NBR dry season 2016');

// ---- Auto-find pixel with the BIGGEST pre/post dieback drop ----
var searchZone = ee.Geometry.Rectangle([135.55, -15.30, 135.90, -14.95]);

var pre = ee.Image(annual.filter(ee.Filter.eq('year', 2014)).first())
  .updateMask(gmw).rename('NBR_pre');
var post = ee.Image(annual.filter(ee.Filter.eq('year', 2016)).first())
  .updateMask(gmw).rename('NBR_post');

var change = pre.subtract(post).rename('drop');   // big positive = big loss
var stack = pre.addBands(post).addBands(change);

// Only consider pixels that were genuinely healthy before (pre > 0.5)
var healthyBefore = stack.updateMask(pre.gt(0.5));

var samples = healthyBefore.sample({
  region: searchZone,
  scale: 30,
  numPixels: 3000,
  geometries: true
});

var sorted = samples.sort('drop', false);   // biggest drop first
var worstDrop = sorted.first();
print('Biggest pre->post drop (feature):', worstDrop);

var diebackPoint = worstDrop.geometry();
Map.centerObject(diebackPoint, 14);
Map.addLayer(diebackPoint, {color: 'red'}, 'Auto-found dieback pixel v2');

var chart = ui.Chart.image.series({
  imageCollection: annual.select('NBR'),
  region: diebackPoint,
  reducer: ee.Reducer.first(),
  scale: 30
}).setOptions({
  title: 'NBR 1987–2025 — biggest 2014-to-2016 drop',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'},
  lineWidth: 2, pointSize: 4
});
print(chart);