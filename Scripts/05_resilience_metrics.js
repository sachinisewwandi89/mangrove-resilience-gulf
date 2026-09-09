// 05_resilience_metrics  —  PHASE 3

// Four resilience surfaces for the 2015-16 dieback, full coast:
//   magnitude of loss · recovery rate · recovery ratio · years to recovery
// plus an unrecovered flag, ENSO-adjusted variants, and aggregation to
// patch and coastal-segment scale.

// =====================================================================
// PART 1 — PHASE 1 PIPELINE  
// =====================================================================
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);
Map.centerObject(gulf, 7);
 
var gmwUnion = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/union/gmw_v3_mng_union');
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
 
var Y0 = 1987, Y1 = 2025;
var years = ee.List.sequence(Y0, Y1);
 
var annual = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    var sub   = allNBR.filterDate(start, end);
    return sub.median().rename('NBR')
      .updateMask(mangroveMask)
      .clip(gulf)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);
print('Composites (should be 39):', annual.size());
 
 
// =====================================================================
// PART 2 — SEGMENTATION  
// =====================================================================
var ltInput = annual.select('NBR').map(function (img) {
  return img.multiply(-1).rename('NBR')
            .set('system:time_start', img.get('system:time_start'));
});
 
var lt = ee.Algorithms.TemporalSegmentation.LandTrendr({
  timeSeries:             ltInput,
  maxSegments:            6,
  spikeThreshold:         0.9,
  vertexCountOvershoot:   3,
  preventOneYearRecovery: true,
  recoveryThreshold:      0.25,
  pvalThreshold:          0.05,
  bestModelProportion:    0.75,
  minObservationsNeeded:  6
});
 
 
// =====================================================================
// PART 3 — FITTED TRAJECTORY, ONE BAND PER YEAR
// =====================================================================

var ASSET_DIR = 'projects/mangroveresilience2026/assets/mangrove-resilience/';
 
// 30 m grid 
var BASE30 = ee.Projection('EPSG:4326').atScale(30);
 
// ---------------------------------------------------------------------
// STAGE SWITCH 
// ---------------------------------------------------------------------

var USE_FITTED_ASSET = true;
 
var yearStrs = [];
for (var y = Y0; y <= Y1; y++) { yearStrs.push('yr_' + y); }
 
var fitted, rmse;
if (USE_FITTED_ASSET) {
  var fa = ee.Image(ASSET_DIR + 'fitted_trajectories');
  fitted = fa.select(yearStrs).updateMask(mangroveMask);
  rmse   = fa.select('rmse').updateMask(mangroveMask);
} else {
  var ltArr  = lt.select('LandTrendr');
  var yrRow  = ltArr.arraySlice(0, 0, 1).arrayProject([1]);
  var fitRow = ltArr.arraySlice(0, 2, 3).arrayProject([1]);
 
  var fittedImgs = [];
  for (var k = 0; k < yearStrs.length; k++) {
    var v = fitRow.arrayMask(yrRow.eq(Y0 + k))
                  .arrayPad([1], -9999)
                  .arraySlice(0, 0, 1)
                  .arrayFlatten([['v']]);
    fittedImgs.push(v.updateMask(v.neq(-9999))
                     .multiply(-1)          // back to NBR orientation
                     .rename(yearStrs[k]));
  }
  fitted = ee.Image.cat(fittedImgs).updateMask(mangroveMask);
  rmse = lt.select('rmse').updateMask(mangroveMask);
}
 
 
// =====================================================================
// PART 4 — TARGET DISTURBANCE
// =====================================================================
var classImg   = ee.Image(ASSET_DIR + 'classMap_fullcoast').select('class');
var diebackImg = ee.Image(ASSET_DIR + 'dieback_2015_16');
 
// Metrics are computed for the 2015-16 dieback (class 2). 
var target = classImg.eq(2).selfMask();
var yod    = diebackImg.select('yod').updateMask(target);
 

function fittedAt(yrImage) {
  // yrImage: per-pixel year
  var stack = ee.ImageCollection(yearStrs.map(function (b, i) {
    var yr = Y0 + i;
    return fitted.select(b).updateMask(yrImage.eq(yr)).rename('v');
  }));
  return stack.mosaic();
}
 
var preNBR = fittedAt(yod.subtract(1)).rename('pre_nbr');
 
 
// =====================================================================
// PART 5 — POST-DISTURBANCE TROUGH
// =====================================================================
 
var troughStack = ee.ImageCollection(yearStrs.map(function (b, i) {
  var yr = Y0 + i;
  return fitted.select(b)
    .updateMask(ee.Image.constant(yr).gte(yod).and(
                ee.Image.constant(yr).lte(yod.add(2))))
    .rename('v');
}));
 
var minNBR = troughStack.min().rename('min_nbr').updateMask(target);
 
var yrMinStack = ee.ImageCollection(yearStrs.map(function (b, i) {
  var yr = Y0 + i;
  return fitted.select(b)
    .updateMask(ee.Image.constant(yr).gte(yod)
        .and(ee.Image.constant(yr).lte(yod.add(2))))
    .eq(minNBR)
    .multiply(yr).toInt16().selfMask().rename('v');
}));
var yearOfMin = yrMinStack.min().rename('year_min').updateMask(target);
 
 
// =====================================================================
// PART 6 — ENSO / REFERENCE-STAND NORMALISATION  
// =====================================================================
 
var REF_SCALE = 5000;
 
var refNBR = fitted.updateMask(classImg.eq(1))
  .setDefaultProjection(BASE30)
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535,
                     bestEffort: true})
  .reproject({crs: 'EPSG:4326', scale: REF_SCALE});
 
 
// Difference-in-differences: what the reference stands did over the same
// interval is subtracted from what the disturbed pixel did. Reference
// values are taken at the same two dates as the pixel's own comparison.
var refAtPre = ee.ImageCollection(yearStrs.map(function (b, i) {
  var yr = Y0 + i;
  return refNBR.select(b).updateMask(yod.subtract(1).eq(yr)).rename('v');
})).mosaic().rename('ref_pre');
 
var refAtNow = refNBR.select('yr_2025').rename('ref_now');
 
var climateShift = refAtNow.subtract(refAtPre).rename('climate_shift');
 
// Coverage of the adjustment
var hasRef = refAtPre.mask().and(refAtNow.mask()).rename('has_ref');
 
 
// =====================================================================
// PART 7 — THE FOUR RESILIENCE MEASURES
// =====================================================================

 
var nowNBR = fitted.select(['yr_2023', 'yr_2024', 'yr_2025'])
  .reduce(ee.Reducer.mean()).rename('now_nbr').updateMask(target);
 
// --- 1. MAGNITUDE OF LOSS -------------------------------------------
var magnitude = preNBR.subtract(minNBR).rename('magnitude');
 
// --- 2. RECOVERY RATIO ----------------------------------------------

var recoveryRatio = nowNBR.subtract(minNBR)
  .divide(magnitude)
  .rename('recovery_ratio');
 
// Climate-adjusted: the reference stands' own change removed.
var recoveryRatioAdj = nowNBR.subtract(climateShift).subtract(minNBR)
  .divide(magnitude)
  .rename('recovery_ratio_adj');
 
// --- 3. RECOVERY RATE -----------------------------------------------
// NBR units regained per year since the trough.
var elapsed = ee.Image.constant(2024).subtract(yearOfMin).max(1);
var recoveryRate = nowNBR.subtract(minNBR).divide(elapsed)
  .rename('recovery_rate');
 
var recoveryRateAdj = nowNBR.subtract(climateShift).subtract(minNBR)
  .divide(elapsed).rename('recovery_rate_adj');
 
// --- 4. YEARS TO RECOVERY -------------------------------------------

var RECOVERY_FRAC = 0.80;
 
var recTarget = minNBR.add(magnitude.multiply(RECOVERY_FRAC));
 
var recYearStack = ee.ImageCollection(yearStrs.map(function (b, i) {
  var yr = Y0 + i;
  return fitted.select(b).gte(recTarget)
    .and(ee.Image.constant(yr).gte(yearOfMin))
    .multiply(yr).toInt16().selfMask().rename('v');
}));
 
var yearRecovered = recYearStack.min().rename('year_recovered')
  .updateMask(target);
 
var yearsToRecovery = yearRecovered.subtract(yod)
  .rename('years_to_recovery');
 
// --- UNRECOVERED FLAG -----------------------------------------------
// 1 where the threshold is never reached by 2025.
var unrecovered = yearRecovered.mask().not()
  .updateMask(target).rename('unrecovered');
 
 
// =====================================================================
// PART 8 — THE RESILIENCE STACK
// =====================================================================
var resilience = magnitude
  .addBands(recoveryRatio)
  .addBands(recoveryRatioAdj)
  .addBands(recoveryRate)
  .addBands(recoveryRateAdj)
  .addBands(yearsToRecovery)
  .addBands(unrecovered)
  .addBands(preNBR)
  .addBands(minNBR)
  .addBands(nowNBR)
  .addBands(yod.rename('yod'))
  .addBands(yearOfMin)
  .addBands(climateShift)
  .addBands(hasRef)
  .addBands(rmse.rename('rmse'))
  .updateMask(target)
  .toFloat();
 
print('Resilience stack bands:', resilience.bandNames());
 
 
// =====================================================================
// PART 9 — AGGREGATION TO PATCH AND COASTAL SEGMENT
// =====================================================================

 
// PATCH: connected components of the dieback class.
var patchId = target.connectedComponents({
  connectedness: ee.Kernel.plus(1),
  maxSize: 1024
}).select('labels').rename('patch_id');
 
var patchMean = resilience.select(
    ['magnitude', 'recovery_ratio', 'recovery_ratio_adj',
     'recovery_rate', 'years_to_recovery', 'unrecovered'])
  .addBands(patchId)
  .reduceConnectedComponents({
    reducer: ee.Reducer.mean(),
    labelBand: 'patch_id',
    maxSize: 1024
  });
 
// COASTAL SEGMENT
var segment = resilience.select(
    ['magnitude', 'recovery_ratio', 'recovery_ratio_adj',
     'recovery_rate', 'years_to_recovery', 'unrecovered'])
  .setDefaultProjection(BASE30)
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535,
                     bestEffort: true})
  .reproject({crs: 'EPSG:4326', scale: REF_SCALE});
 
 
// =====================================================================
// PART 10 — DIAGNOSTICS 
// =====================================================================
var STATS = ee.Reducer.mean()
  .combine(ee.Reducer.median(), '', true)
  .combine(ee.Reducer.stdDev(), '', true)
  .combine(ee.Reducer.percentile([10, 90]), '', true);
 
print('R1 — dieback area (ha):',
  ee.Image.pixelArea().divide(10000).updateMask(target).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: gulf, scale: 30,
    maxPixels: 1e13, tileScale: 4}));
 
print('R2 — metric distributions over the dieback:',
  resilience.select(['magnitude', 'recovery_ratio', 'recovery_ratio_adj',
                     'recovery_rate', 'years_to_recovery'])
    .reduceRegion({reducer: STATS, geometry: gulf, scale: 90,
                   maxPixels: 1e13, tileScale: 16}));
 
// THE HEADLINE NUMBER: what fraction of the 2015-16 dieback has not
// recovered to 80% of pre-disturbance condition after nine years.
print('R3 — unrecovered fraction (mean of the 0/1 flag):',
  unrecovered.reduceRegion({reducer: ee.Reducer.mean(), geometry: gulf,
    scale: 30, maxPixels: 1e13, tileScale: 8}));
 
// THE SIZE OF THE ENSO CORRECTION. If climate_shift is large relative
// to magnitude, the raw recovery ratio is mostly climate and decision
// 15 is vindicated. If it is small, say so — a null result on a
// documented decision is still a result.
print('R4 — climate shift vs magnitude:',
  climateShift.rename('shift').addBands(magnitude)
    .reduceRegion({reducer: STATS, geometry: gulf, scale: 90,
                   maxPixels: 1e13, tileScale: 16}));
 
print('R5 — reference coverage (fraction of dieback with a 5 km reference):',
  hasRef.updateMask(target).unmask(0).updateMask(target)
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: gulf, scale: 90,
                   maxPixels: 1e13, tileScale: 8}));
 
Map.addLayer(resilience.select('recovery_ratio_adj'),
  {min: 0, max: 1.2, palette: ['67001f', 'f7f7f7', '1a9850']},
  'Recovery ratio (ENSO-adjusted)', true);
Map.addLayer(resilience.select('unrecovered'),
  {min: 0, max: 1, palette: ['ffffff', 'd73027']}, 'Unrecovered', false);
Map.addLayer(resilience.select('years_to_recovery'),
  {min: 1, max: 9, palette: ['1a9850', 'fee08b', 'd73027']},
  'Years to recovery', false);
 
 
// =====================================================================
// PART 11 — EXPORTS   (start these as early in the day as possible)
// =====================================================================
Export.image.toAsset({
  image: resilience,
  description: 'resilience_2015_16',
  assetId: ASSET_DIR + 'resilience_2015_16',
  region: gulf, scale: 30, maxPixels: 1e13, crs: 'EPSG:4326'
});
 
Export.image.toAsset({
  image: patchMean,
  description: 'resilience_patch',
  assetId: ASSET_DIR + 'resilience_patch',
  region: gulf, scale: 30, maxPixels: 1e13, crs: 'EPSG:4326'
});
 
Export.image.toAsset({
  image: segment,
  description: 'resilience_segment',
  assetId: ASSET_DIR + 'resilience_segment',
  region: gulf, scale: REF_SCALE, maxPixels: 1e13, crs: 'EPSG:4326'
});
 
// STAGE A TASK — start this one FIRST and alone. 40 bands: the 39 fitted
// years plus rmse. Once it exists, set USE_FITTED_ASSET = true and every
// other task above reads it instead of re-running the segmentation.
// Phase 4 reads it too, so the segmentation is computed exactly once for
// the whole project.
if (!USE_FITTED_ASSET) {
  Export.image.toAsset({
    image: fitted.addBands(rmse.rename('rmse')).toFloat(),
    description: 'fitted_trajectories',
    assetId: ASSET_DIR + 'fitted_trajectories',
    region: gulf, scale: 30, maxPixels: 1e13, crs: 'EPSG:4326'
  });
}
 
 
// PART 12 — RESULTS PACKAGE
//
// The prints in Part 10 will very likely time out: reduceRegion over
// 920,000 pixels in the interactive console has a 5-minute ceiling.
// Everything here runs as a BATCH task instead, which has no such
// limit, and lands as CSV in Drive ready for figures.
//
// Four outputs:
//   phase3_trajectories.csv  39 rows — the headline figure
//   phase3_summary.csv       1 row  — every number for the results text
//   phase3_pixels.csv        <=50k rows — distributions and scatter plots
//   phase3_segments.csv      5 km coastal segments — the map figure
//
// Start all four now. They are far cheaper than the image exports.
// =====================================================================
 
 
// ---------------------------------------------------------------------
// 12.1 — MEAN TRAJECTORIES, DISTURBED vs UNDISTURBED
// ---------------------------------------------------------------------
// This is Figure 1 of the manuscript. Two lines, 1987-2025: the mean
// fitted NBR of the 2015-16 dieback pixels, and the mean of undisturbed
// mangrove in the same region. The gap between them after 2016 is the
// unrecovered deficit, and the fact that BOTH lines rise into 2021-22
// is the visual proof that the ENSO correction is necessary.
//
// scale 90 rather than 30: this is a regional mean over hundreds of
// thousands of pixels and the third decimal place will not move.
 
var refClass1 = classImg.eq(1);
 
var trajectories = ee.FeatureCollection(yearStrs.map(function (b, i) {
  var yr = Y0 + i;
  var rawYr = ee.Image(annual.filter(ee.Filter.eq('year', yr)).first())
                .select('NBR');
 
  var fitDis = fitted.select(b).updateMask(target).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: gulf, scale: 90,
    maxPixels: 1e13, tileScale: 16});
  var fitRef = fitted.select(b).updateMask(refClass1).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: gulf, scale: 90,
    maxPixels: 1e13, tileScale: 16});
  var rawDis = rawYr.updateMask(target).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: gulf, scale: 90,
    maxPixels: 1e13, tileScale: 16});
  var rawRef = rawYr.updateMask(refClass1).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: gulf, scale: 90,
    maxPixels: 1e13, tileScale: 16});
 
  return ee.Feature(null, {
    year:          yr,
    fit_dieback:   fitDis.get(b),
    fit_reference: fitRef.get(b),
    raw_dieback:   rawDis.get('NBR'),
    raw_reference: rawRef.get('NBR')
  });
}));
 
Export.table.toDrive({
  collection: trajectories,
  description: 'phase3_trajectories',
  folder: 'mangrove_resilience',
  fileFormat: 'CSV',
  selectors: ['year', 'fit_dieback', 'fit_reference',
              'raw_dieback', 'raw_reference']
});
 
 
// ---------------------------------------------------------------------
// 12.2 — SUMMARY ROW
// ---------------------------------------------------------------------
// Every number that goes in the results paragraph, in one row.
//
// NOTE ON years_to_recovery: it is masked where recovery never happens,
// so its mean is the mean AMONG PIXELS THAT RECOVERED. That is the
// correct quantity, but it must be reported that way — "of the pixels
// that recovered, median X years" — never as a mean over the dieback.
 
var SUMSTATS = ee.Reducer.mean()
  .combine(ee.Reducer.median(), '', true)
  .combine(ee.Reducer.stdDev(), '', true)
  .combine(ee.Reducer.percentile([10, 25, 75, 90]), '', true);
 
var METRICS = ['magnitude', 'recovery_ratio', 'recovery_ratio_adj',
               'recovery_rate', 'recovery_rate_adj', 'years_to_recovery',
               'pre_nbr', 'min_nbr', 'now_nbr', 'climate_shift', 'rmse'];
 
var sumArea = ee.Image.pixelArea().divide(10000).updateMask(target)
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: gulf, scale: 30,
                 maxPixels: 1e13, tileScale: 8});
 
var sumPix = target.reduceRegion({
  reducer: ee.Reducer.count(), geometry: gulf, scale: 30,
  maxPixels: 1e13, tileScale: 8});
 
var sumMetrics = resilience.select(METRICS).reduceRegion({
  reducer: SUMSTATS, geometry: gulf, scale: 30,
  maxPixels: 1e13, tileScale: 16});
 
var sumUnrec = unrecovered.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: gulf, scale: 30,
  maxPixels: 1e13, tileScale: 8});
 
// Unrecovered area in hectares — the number a manager acts on.
var sumUnrecHa = ee.Image.pixelArea().divide(10000)
  .updateMask(unrecovered.eq(1))
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: gulf, scale: 30,
                 maxPixels: 1e13, tileScale: 8});
 
var sumRefCov = hasRef.unmask(0).updateMask(target).reduceRegion({
  reducer: ee.Reducer.mean(), geometry: gulf, scale: 90,
  maxPixels: 1e13, tileScale: 8});
 
// Sensitivity of the headline to the 80% threshold. Three thresholds,
// so the results section can state that the unrecovered fraction is
// robust (or is not) rather than resting on one arbitrary choice.
function unrecAt(frac) {
  var tgt = minNBR.add(magnitude.multiply(frac));
  var got = ee.ImageCollection(yearStrs.map(function (b, i) {
    var yr = Y0 + i;
    return fitted.select(b).gte(tgt)
      .and(ee.Image.constant(yr).gte(yearOfMin))
      .multiply(yr).toInt16().selfMask().rename('v');
  })).min();
  return got.mask().not().updateMask(target).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: gulf, scale: 30,
    maxPixels: 1e13, tileScale: 8}).values().get(0);
}
 
var summary = ee.Feature(null, ee.Dictionary({
    dieback_ha:        sumArea.get('area'),
    dieback_pixels:    sumPix.get('class'),
    unrecovered_frac:  sumUnrec.get('unrecovered'),
    unrecovered_ha:    sumUnrecHa.get('area'),
    ref_coverage:      sumRefCov.get('has_ref'),
    unrec_frac_70:     unrecAt(0.70),
    unrec_frac_80:     unrecAt(0.80),
    unrec_frac_90:     unrecAt(0.90)
  }).combine(sumMetrics));
 
Export.table.toDrive({
  collection: ee.FeatureCollection([summary]),
  description: 'phase3_summary',
  folder: 'mangrove_resilience',
  fileFormat: 'CSV'
});
 
 
// ---------------------------------------------------------------------
// 12.3 — PER-PIXEL SAMPLE
// ---------------------------------------------------------------------
// For histograms, the magnitude-versus-recovery scatter, and anything
// else that needs a distribution rather than a summary. Do the plotting
// in Python — GEE charting will not produce publication figures.
//
// CRITICAL: unmask BEFORE sampling. years_to_recovery is masked exactly
// where recovery never happened, so sampling with dropNulls would
// silently delete every unrecovered pixel — the most important half of
// the result. -9999 is the sentinel; filter on it in Python.
 
var forSample = resilience.unmask(-9999).updateMask(target)
  .addBands(ee.Image.pixelLonLat());
 
var pixelSample = forSample.sample({
  region:    gulf,
  scale:     30,
  numPixels: 50000,
  seed:      2026,
  tileScale: 16,
  dropNulls: false,
  geometries: false
});
 
Export.table.toDrive({
  collection: pixelSample,
  description: 'phase3_pixels',
  folder: 'mangrove_resilience',
  fileFormat: 'CSV'
});
 
 
// ---------------------------------------------------------------------
// 12.4 — COASTAL SEGMENTS
// ---------------------------------------------------------------------
// 5 km segment means with coordinates: the spatial figure, and the unit
// at which Phase 4 blocking will operate. Exporting it now means the
// blocking scheme is already defined when Phase 4 starts.
 
var segPts = segment.addBands(ee.Image.pixelLonLat())
  .addBands(ee.Image.pixelArea().divide(10000).rename('cell_ha'))
  .sample({
    region:     gulf,
    scale:      REF_SCALE,
    projection: 'EPSG:4326',
    numPixels:  20000,
    seed:       2026,
    tileScale:  8,
    dropNulls:  true,
    geometries: false
  });
 
Export.table.toDrive({
  collection: segPts,
  description: 'phase3_segments',
  folder: 'mangrove_resilience',
  fileFormat: 'CSV'
});