// 03_change_detection  —  investigating 2011, and the D7 mask test

// THE QUESTION THIS ANSWERS
// LandTrendr found 3,619 detections in 2011 - 35% of everything, and
// the modal detection year. But 2010 was the highest NBR year in the
// entire 39-year record (record La Nina). So: did something actually
// happen in 2011, or did the forest simply return to normal after an
// exceptionally wet year?

// =====================================================================
// PART 1 — PHASE 1 PIPELINE 
// =====================================================================
 
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);

Map.centerObject(gulf, 7);
 
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
  var nbr = sr.normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR');
  return ee.Image(nbr).copyProperties(img, ['system:time_start']);
}
 
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
 
var allNBR = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepTMETM)
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepTMETM))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepOLI))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
      .filterBounds(gulf).map(maskL2).map(prepOLI));
 
var years = ee.List.sequence(1987, 2025);
 
var annualUnmasked = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    var sub   = allNBR.filterDate(start, end);
    return sub.median().rename('NBR')
      .addBands(sub.count().rename('n_obs').unmask(0))
      .clip(gulf)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);
 
var annual = annualUnmasked.map(function (img) {
  return img.updateMask(mangroveMask);
});
 
print('Composites (should be 39):', annual.size());
 

var testStrip = ee.Geometry.Rectangle([135.55, -15.30, 135.90, -14.95]);
var aoi = testStrip;

// =====================================================================
// PART 2 — CHANGE DETECTION: WHAT HAPPENED IN 2011?
// =====================================================================
 
// mean NBR over a span of years
function meanNBR(y1, y2) {
  return annual.select('NBR')
    .filter(ee.Filter.and(ee.Filter.gte('year', y1),
                          ee.Filter.lte('year', y2)))
    .mean();
}
 
// ---- Reference periods ----------------------------------------------
var base_2007_09 = meanNBR(2007, 2009);   // stable pre-La Nina baseline
var yr_2010      = meanNBR(2010, 2010);   // the anomalous peak
var yr_2011      = meanNBR(2011, 2011);
var base_2012_14 = meanNBR(2012, 2014);   // pre-dieback baseline
var yr_2016      = meanNBR(2016, 2016);   // post-dieback
 
 
// ---- TEST 1. Was 2010 genuinely anomalous? --------------------------
// If 2010 sits well above the 2007-09 baseline, then any decline in
// 2011 is at least partly a return to normal rather than a loss.
var anomaly2010 = yr_2010.subtract(base_2007_09).rename('anom2010');
 
print('TEST 1 - mean 2010 anomaly vs 2007-09 baseline:',
  anomaly2010.reduceRegion({reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
Map.addLayer(anomaly2010.updateMask(mangroveMask),
  {min: -0.1, max: 0.1, palette: ['0000ff','ffffff','00ff00']},
  '2010 anomaly (green = higher than baseline)', false);
 
 
// ---- TEST 2. ---------------------------------------
// Two differencing operations:
//
//   dNBR_from_peak     = 2010 - 2011   (what LandTrendr responds to)
//   dNBR_from_baseline = 2007-09 - 2011 (whether real loss occurred)

var dNBR_peak = yr_2010.subtract(yr_2011).rename('dNBR_peak');
var dNBR_base = base_2007_09.subtract(yr_2011).rename('dNBR_base');
 
print('TEST 2a - mean drop from 2010 peak to 2011:',
  dNBR_peak.reduceRegion({reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
print('TEST 2b - mean drop from 2007-09 baseline to 2011:',
  dNBR_base.reduceRegion({reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
// How many pixels show REAL loss in 2011 (below the stable baseline)?
var realLoss2011 = dNBR_base.gt(0.15).updateMask(mangroveMask).selfMask();
print('TEST 2c - pixels below 2007-09 baseline by >0.15 in 2011:',
  realLoss2011.reduceRegion({reducer: ee.Reducer.count(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
Map.addLayer(realLoss2011, {palette: ['ff00ff']},
  '2011: real loss vs baseline', false);
 
 
// ---- TEST 3. The same test for 2015-16, as a positive control -------

var dNBR_dieback = base_2012_14.subtract(yr_2016).rename('dNBR_dieback');
 
print('TEST 3a - mean drop from 2012-14 baseline to 2016:',
  dNBR_dieback.reduceRegion({reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
var realLoss2016 = dNBR_dieback.gt(0.15).updateMask(mangroveMask).selfMask();
print('TEST 3b - pixels below 2012-14 baseline by >0.15 in 2016:',
  realLoss2016.reduceRegion({reducer: ee.Reducer.count(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
Map.addLayer(realLoss2016, {palette: ['ff0000']},
  '2016: real loss vs baseline', true);
 
 
// ---- TEST 4. Did 2011 recover? --------------------------------------

var recovery2011 = base_2012_14.subtract(base_2007_09).rename('rec2011');
 
print('TEST 4 - 2012-14 mean minus 2007-09 mean (near zero = recovered):',
  recovery2011.reduceRegion({reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 30, maxPixels: 1e13}));
 
// =====================================================================
// PART 3 — D7: PRE-1996 MASK SURVIVORSHIP TEST
// =====================================================================
// The GMW mask covers 1996-2020. A stand that existed in 1987, died
// before 1996 and never returned is absent from every epoch, therefore
// absent from the mask, therefore invisible to LandTrendr.
//
// The bias has a direction: the pre-1996 disturbance record is skewed
// toward events the mangrove RECOVERED from. 
//
// This matters more since Duke et al. (2022) identified a 1982 dieback
// with ~15-year recovery, which would place regrowth right across the
// early part of this series.
 
var unionBin = gmwUnion.unmask(0).gt(0);
 
var buffer = unionBin.focalMax(30, 'square', 'pixels')
                     .and(unionBin.not())
                     .selfMask().clip(gulf);
 
var early = annualUnmasked.select('NBR')
  .filter(ee.Filter.and(ee.Filter.gte('year', 1987),
                        ee.Filter.lte('year', 1992))).median();
 
var late = annualUnmasked.select('NBR')
  .filter(ee.Filter.and(ee.Filter.gte('year', 2015),
                        ee.Filter.lte('year', 2020))).median();
 
// Candidates
var candidates = early.gt(0.45).and(late.lt(0.25))
                      .updateMask(buffer).selfMask();
 
print('D7 - candidate lost-mangrove area (ha):',
  ee.Image.pixelArea().divide(10000).updateMask(candidates)
    .reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi,
                   scale: 30, maxPixels: 1e13}));
 
Map.addLayer(candidates, {palette: ['ff00ff']},
  'D7: candidate pre-1996 losses', false);
 
