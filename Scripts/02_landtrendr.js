// 02_landtrendr  —  PHASE 2: disturbance detection + validation setup

// PART 1 — PHASE 1 PIPELINE 
 
// ---- 1. Study area --------------------------------------------------
// Southern and western Gulf, NT + QLD
// Northern limit 14.0 S: north of this the coast becomes northeast
// Arnhem Land, which is outside the dieback region.
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);
 
var testStrip = ee.Geometry.Rectangle([135.55, -15.30, 135.90, -14.95]);
var aoi = testStrip;
 
 
// ---- 2. Mangrove mask -----------------------------------------------
var gmwUnion = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/union/gmw_v3_mng_union');
var gmwCore = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/core/gmw_v3_mng_core');
 
var mangroveMask = gmwUnion.gt(0).selfMask().clip(gulf);
 
 
// ---- 3. Cloud / shadow / fill mask 
function maskL2(img) {
  var qa = img.select('QA_PIXEL');
  var bad = qa.bitwiseAnd(1 << 0).neq(0)      // fill
      .or(qa.bitwiseAnd(1 << 1).neq(0))       // dilated cloud
      .or(qa.bitwiseAnd(1 << 2).neq(0))       // cirrus
      .or(qa.bitwiseAnd(1 << 3).neq(0))       // cloud
      .or(qa.bitwiseAnd(1 << 4).neq(0));      // cloud shadow
  return img.updateMask(bad.not());
}
 
 
// ---- 4. Roy et al. (2016) harmonisation

var ROY = {
  nir:   {slope: 0.8462, itcp: 0.0412},
  swir2: {slope: 0.9071, itcp: 0.0172}
};
 
function scaleSR(img) {
  return img.select('SR_B.').multiply(0.0000275).add(-0.2);
}
 
// L5 / L7: NIR = SR_B4, SWIR2 = SR_B7
function prepTMETM(img) {
  var sr = scaleSR(img);
  var nbr = sr.normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR');
  return ee.Image(nbr).copyProperties(img, ['system:time_start']);
}
 
// L8 / L9: NIR = SR_B5, SWIR2 = SR_B7. Harmonise to ETM+ first.
function prepOLI(img) {
  var sr = scaleSR(img);
  var nir = sr.select('SR_B5')
              .multiply(ROY.nir.slope).add(ROY.nir.itcp).rename('NIR');
  var swir2 = sr.select('SR_B7')
              .multiply(ROY.swir2.slope).add(ROY.swir2.itcp).rename('SWIR2');
  var nbr = nir.addBands(swir2)
               .normalizedDifference(['NIR', 'SWIR2']).rename('NBR');
  return ee.Image(nbr).copyProperties(img, ['system:time_start']);
}
 
 
// ---- 5. Merged NBR collection 
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepOLI);
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepOLI);
 
var allNBR = l5.merge(l7).merge(l8).merge(l9);
 
 
// ---- 6. Annual dry-season composites 

var years = ee.List.sequence(1987, 2025);
 
var annual = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    var sub   = allNBR.filterDate(start, end);
 
    var med  = sub.median().rename('NBR');
    var nObs = sub.count().rename('n_obs').unmask(0);
 
    return med.addBands(nObs)
      .updateMask(mangroveMask)
      .clip(gulf)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);
 
print('Composites :', annual.size());


// PART 2 — LANDTRENDR SEGMENTATION

 
// ---- 7. Run LandTrendr on harmonised composites 

var ltInput = annual.select('NBR').map(function (img) {
  return img.multiply(-1).rename('NBR')
            .set('system:time_start', img.get('system:time_start'));
});
 
var lt = ee.Algorithms.TemporalSegmentation.LandTrendr({
  timeSeries:             ltInput.filterBounds(aoi),
  maxSegments:            6,
  spikeThreshold:         0.9,
  vertexCountOvershoot:   3,
  preventOneYearRecovery: true,
  recoveryThreshold:      0.25,
  pvalThreshold:          0.05,
  bestModelProportion:    0.75,
  minObservationsNeeded:        6
});
 
print('LandTrendr output (expect LandTrendr + rmse bands):', lt);


// ---- 8. Extract the greatest disturbance per pixel 
var ltgee = require('users/emaprlab/public:Modules/LandTrendr.js');
 
var changeParams = {
  delta:  'loss',
  sort:   'greatest',
  year:   {checked: true,  start: 1987, end: 2025},
  mag:    {checked: true,  value: 0.15, operator: '>'},
  dur:    {checked: false, value: 4,   operator: '<'},
  preval: {checked: false, value: 300, operator: '>'},
  mmu:    {checked: true,  value: 5}
};

// The module uses this to set the flip direction when it
// converts fitted magnitudes back to the original index orientation.
changeParams.index = 'NBR';

var changeImg = ltgee.getChangeMap(lt, changeParams).clip(aoi);

var yod = changeImg.select('yod').updateMask(mangroveMask);
var mag = changeImg.select('mag').updateMask(mangroveMask);

// ---- DIAGNOSTIC A: is the mangrove mask attaching? 
print('Mangrove pixels in strip:', mangroveMask.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 30, maxPixels: 1e13}));

// ---- DIAGNOSTIC B: full year distribution of detections 
print('Detections per year:', yod.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(),
  geometry: aoi, scale: 30, maxPixels: 1e13}));

// ---- DIAGNOSTIC C: greatest loss WITHIN 2015-16 only ----
 
var dieback = ltgee.getChangeMap(lt, {
  delta:  'loss',
  sort:   'greatest',
  year:   {checked: true,  start: 2015, end: 2016},
  mag:    {checked: true,  value: 0.15, operator: '>'},
  dur:    {checked: false, value: 4,   operator: '<'},
  preval: {checked: false, value: 300, operator: '>'},
  mmu:    {checked: true,  value: 5},
  index:  'NBR'
}).clip(aoi);

print('2015-16 detections:', dieback.select('yod').reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 30, maxPixels: 1e13}));

Map.addLayer(dieback.select('yod'), {min: 2015, max: 2016,
  palette: ['ff0000','ffaa00']}, '2015-16 dieback only');

print('changeImg band names:', changeImg.bandNames());

// Is anything detected BEFORE masking to mangrove?
print('Detections before mangrove mask:', changeImg.select('yod')
  .reduceRegion({reducer: ee.Reducer.count(), geometry: aoi,
                 scale: 30, maxPixels: 1e13}));

// And after?
print('Detections after mangrove mask:', yod
  .reduceRegion({reducer: ee.Reducer.count(), geometry: aoi,
                 scale: 30, maxPixels: 1e13}));

// Range of detection years, ignoring the mask
print('Year range (unmasked):', changeImg.select('yod')
  .reduceRegion({reducer: ee.Reducer.minMax(), geometry: aoi,
                 scale: 30, maxPixels: 1e13}));

// Bands: yod (year of detection), mag, dur, preval, rate, dsnr
Map.centerObject(aoi, 11);
Map.addLayer(yod, {min: 1990, max: 2025,
  palette: ['blue','green','yellow','orange','red']}, 'Year of detection');
Map.addLayer(mag, {min: 0, max: 800,
  palette: ['ffffff','ff0000']}, 'Magnitude of loss');
 
// Does the dieback date to 2015->2016?

print('Modal detection year in strip:', yod.reduceRegion({
  reducer: ee.Reducer.mode(), geometry: aoi, scale: 30, maxPixels: 1e13
}));

// ---- 9. THE THREE-CLASS MAP ----------------------------------------
// 1 = no disturbance detected 1987-2025
// 2 = disturbance detected 2015-16  (from the WINDOWED run)
// 3 = disturbance detected in some other year

var diebackYod = dieback.select('yod').updateMask(mangroveMask);

var anyDetection = yod.gt(0).unmask(0);
var inWindow     = diebackYod.gt(0).unmask(0);

var classMap = ee.Image(1)
  .where(anyDetection.eq(1), 3)          // any disturbance, other year
  .where(inWindow.eq(1), 2)              // 2015-16 overrides
  .rename('class')
  .updateMask(mangroveMask)
  .clip(aoi)
  .toByte();

Map.addLayer(classMap, {min: 1, max: 3,
  palette: ['cccccc', 'ff0000', 'ffaa00']}, '3-class disturbance map');
  
  // ---- 10. Class areas -> weights for Pontius equation (1) ------------

print('Area (ha) by class:',
  ee.Image.pixelArea().divide(10000).addBands(classMap).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: aoi, scale: 30, maxPixels: 1e13
  }));
  
 // ---- 11. Detection confidence by era 
// Era boundaries follow Landsat sensor availability:
//   1987-1998  L5 only   ~10 obs/season
//   1999-2012  L5 + L7   ~15 obs/season
//   2013-2025  L8 / L9   ~26 obs/season

var era = ee.Image(0)
  .where(yod.gte(1987).and(yod.lte(1998)), 1)
  .where(yod.gte(1999).and(yod.lte(2012)), 2)
  .where(yod.gte(2013).and(yod.lte(2025)), 3)
  .rename('era')
  .updateMask(yod.mask());
 
print('Detected disturbance area (ha) by era:',
  ee.Image.pixelArea().divide(10000).addBands(era).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'era'}),
    geometry: aoi, scale: 30, maxPixels: 1e13
  }));
