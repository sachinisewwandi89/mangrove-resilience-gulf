// 04_full_coast  —  PHASE 2 PRODUCTION RUN
//
// Full-coast LandTrendr segmentation, 920,278 mangrove pixels,
// 1987-2025. 
// =====================================================================
// PART 1 — PHASE 1 PIPELINE
// =====================================================================
 
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);
Map.centerObject(gulf, 7);
 
// FULL COAST
var aoi = gulf;
 
var gmwUnion = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/union/gmw_v3_mng_union');
var gmwCore = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/core/gmw_v3_mng_core');
 
var mangroveMask = gmwUnion.gt(0).selfMask().clip(gulf);
 
function maskL2(img) {
  var qa = img.select('QA_PIXEL');
  var bad = qa.bitwiseAnd(1 << 0).neq(0)
      .or(qa.bitwiseAnd(1 << 1).neq(0))
      .or(qa.bitwiseAnd(1 << 2).neq(0))
      .or(qa.bitwiseAnd(1 << 3).neq(0))
      .or(qa.bitwiseAnd(1 << 4).neq(0));
  return img.updateMask(bad.not());
}
 
var ROY = {
  nir:   {slope: 0.8462, itcp: 0.0412},
  swir2: {slope: 0.9071, itcp: 0.0172}
};
 
function scaleSR(img) {
  return img.select('SR_B.').multiply(0.0000275).add(-0.2);
}
 
function prepTMETM(img) {
  var sr = scaleSR(img);
  return ee.Image(sr.normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR'))
    .copyProperties(img, ['system:time_start']);
}
 
function prepOLI(img) {
  var sr = scaleSR(img);
  var nir = sr.select('SR_B5')
              .multiply(ROY.nir.slope).add(ROY.nir.itcp).rename('NIR');
  var swir2 = sr.select('SR_B7')
              .multiply(ROY.swir2.slope).add(ROY.swir2.itcp).rename('SWIR2');
  return ee.Image(nir.addBands(swir2)
    .normalizedDifference(['NIR', 'SWIR2']).rename('NBR'))
    .copyProperties(img, ['system:time_start']);
}
 
var allNBR = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepTMETM)
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepTMETM))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepOLI))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepOLI));
 
var years = ee.List.sequence(1987, 2025);
 
var annual = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    var sub   = allNBR.filterDate(start, end);
    return sub.median().rename('NBR')
      .addBands(sub.count().rename('n_obs').unmask(0))
      .updateMask(mangroveMask)
      .clip(gulf)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);
 
print('Composites (should be 39):', annual.size());

// =====================================================================
// PART 2 — LANDTRENDR SEGMENTATION  (parameters LOCKED)
// =====================================================================
// NBR is inverted because the algorithm assumes increasing values
// represent disturbance, and NBR falls when canopy is lost.
var ltInput = annual.select('NBR').map(function (img) {
  return img.multiply(-1).rename('NBR')
            .set('system:time_start', img.get('system:time_start'));
});
 
var LT_PARAMS = {
  maxSegments:            6,
  spikeThreshold:         0.9,
  vertexCountOvershoot:   3,
  preventOneYearRecovery: true,
  recoveryThreshold:      0.25,
  pvalThreshold:          0.05,
  bestModelProportion:    0.75,
  minObservationsNeeded:  6
};
 
var lt = ee.Algorithms.TemporalSegmentation.LandTrendr({
  timeSeries:             ltInput,
  maxSegments:            LT_PARAMS.maxSegments,
  spikeThreshold:         LT_PARAMS.spikeThreshold,
  vertexCountOvershoot:   LT_PARAMS.vertexCountOvershoot,
  preventOneYearRecovery: LT_PARAMS.preventOneYearRecovery,
  recoveryThreshold:      LT_PARAMS.recoveryThreshold,
  pvalThreshold:          LT_PARAMS.pvalThreshold,
  bestModelProportion:    LT_PARAMS.bestModelProportion,
  minObservationsNeeded:  LT_PARAMS.minObservationsNeeded
});
 
 
// =====================================================================
// PART 3 — CHANGE EXTRACTION  (thresholds LOCKED)
// =====================================================================
// mag 0.15 is in NATIVE NBR units. It sits above the regional ENSO
// amplitude of ~0.08 measured in D4. Sensitivity tested: 0.15 vs 0.10
// changed dieback area by 2.4%.
//
// mmu 5 pixels = 0.45 ha. Accad et al. (2019) mapped 653 Queensland
// patches averaging 4 ha, smallest 0.02 ha. Undercount of the smallest
// is structural at 30 m.
var ltgee = require('users/emaprlab/public:Modules/LandTrendr.js');
 
var baseParams = {
  delta:  'loss',
  sort:   'greatest',
  mag:    {checked: true,  value: 0.15, operator: '>'},
  dur:    {checked: false, value: 4,    operator: '<'},
  preval: {checked: false, value: 300,  operator: '>'},
  mmu:    {checked: true,  value: 5},
  index:  'NBR'
};
 
// PRODUCT A - unrestricted disturbance history 
var histParams = JSON.parse(JSON.stringify(baseParams));
histParams.year = {checked: true, start: 1987, end: 2025};
var history = ltgee.getChangeMap(lt, histParams).clip(gulf);
 
// PRODUCT B - the 2015-16 dieback map

var dbParams = JSON.parse(JSON.stringify(baseParams));
dbParams.year = {checked: true, start: 2015, end: 2016};
var dieback = ltgee.getChangeMap(lt, dbParams).clip(gulf);
 
 
// =====================================================================
// PART 4 — THE THREE-CLASS MAP
// =====================================================================
// 1 = no disturbance detected 1987-2025
// 2 = disturbance detected 2015-16
// 3 = disturbance detected in another year
var histYod = history.select('yod').updateMask(mangroveMask);
var dbYod   = dieback.select('yod').updateMask(mangroveMask);
 
var classMap = ee.Image(1)
  .where(histYod.gt(0).unmask(0).eq(1), 3)
  .where(dbYod.gt(0).unmask(0).eq(1), 2)   // 2015-16 overrides
  .rename('class')
  .updateMask(mangroveMask)
  .clip(gulf)
  .toByte();
 
 
// =====================================================================
// PART 5 — 2011 EXCLUSION FLAG
// =====================================================================
// Change detection established that 2011 detections are
// regression from the 2010 La Nina maximum, not disturbance:
//   2010 sat +0.105 NBR above the 2007-09 baseline
//   2011 was -0.033 BELOW that peak but still ABOVE baseline
//   361 pixels showed real loss against 3,619 segmentation detections
//   2012-14 was back to baseline (+0.019)

var suspectYear = histYod.eq(2011)
  .or(histYod.eq(2001))
  .or(histYod.eq(2013))     // also the L7->L8 transition year
  .rename('suspect')
  .updateMask(mangroveMask);
 
 
// =====================================================================
// PART 6 — VISUAL CHECK BEFORE EXPORTING
// =====================================================================
 Map.addLayer(classMap, {min: 1, max: 3,
   palette: ['cccccc', 'ff0000', 'ffaa00']}, '3-class map', true);
 Map.addLayer(dbYod, {min: 2015, max: 2016,
   palette: ['ff0000', 'ffaa00']}, '2015-16 dieback', false);
 

print('FULL COAST area (ha) by class - STRATUM WEIGHTS:',
  ee.Image.pixelArea().divide(10000).addBands(classMap).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: gulf, scale: 30, maxPixels: 1e13
  }));
 

 // =====================================================================
// PART 7 — EXPORTS   
// =====================================================================

 
Export.image.toAsset({
  image: classMap,
  description: 'classMap_fullcoast',
  assetId: 'users/sachinidd/mangrove-resilience/classMap_fullcoast',
  region: gulf, scale: 30, maxPixels: 1e13, crs: 'EPSG:4326'
});
 
Export.image.toAsset({
  image: dieback.select(['yod','mag','dur','preval','rate','dsnr']),
  description: 'dieback_2015_16',
  assetId: 'users/sachinidd/mangrove-resilience/dieback_2015_16',
  region: gulf, scale: 30, maxPixels: 1e13, crs: 'EPSG:4326'
});
 
Export.image.toAsset({
  image: history.select(['yod','mag','dur','preval','rate','dsnr'])
                .addBands(suspectYear),
  description: 'history_fullcoast',
  assetId: 'users/sachinidd/mangrove-resilience/history_fullcoast',
  region: gulf, scale: 30, maxPixels: 1e13, crs: 'EPSG:4326'
});
 

