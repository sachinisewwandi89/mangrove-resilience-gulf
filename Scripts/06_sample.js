// 06_sample  —  PHASE 2 VALIDATION: STRATIFIED REFERENCE SAMPLE

// =====================================================================
// 1 — STUDY AREA, MASK, ASSETS
// =====================================================================
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);
Map.centerObject(gulf, 7);
 
var gmwUnion = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/union/gmw_v3_mng_union');
var mangroveMask = gmwUnion.gt(0).selfMask().clip(gulf);
 
var ASSET_DIR = 'projects/mangroveresilience2026/assets/mangrove-resilience/';
 
var classImg   = ee.Image(ASSET_DIR + 'classMap_fullcoast');
var diebackImg = ee.Image(ASSET_DIR + 'dieback_2015_16');
var historyImg = ee.Image(ASSET_DIR + 'history_fullcoast');
 
print('classMap bands (expect 1):', classImg.bandNames());
print('dieback bands (expect 6):', diebackImg.bandNames());
print('history bands (expect 7, incl. suspect):', historyImg.bandNames());
 
 
// =====================================================================
// 2 — ANNUAL COMPOSITES  
// =====================================================================
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
 
var nbrBandNames = [];
var nbrImages = [];
for (var y = Y0; y <= Y1; y++) {
  var nm = 'nbr_' + y;
  nbrBandNames.push(nm);
  nbrImages.push(
    allNBR.filterDate(ee.Date.fromYMD(y, 5, 1), ee.Date.fromYMD(y, 10, 31))
          .median().rename(nm).updateMask(mangroveMask)
  );
}
var nbrStack = ee.Image.cat(nbrImages).unmask(-9999);
 
 
// =====================================================================
// 3 — STRATIFICATION 
// =====================================================================


 
var SPLIT_CLASS3 = false;      
 
var SEED = 2026;
 
var CLASS_VALUES, CLASS_POINTS, STRATUM_NOTE;
if (SPLIT_CLASS3) {
  CLASS_VALUES = [1, 2, 3, 4];
  CLASS_POINTS = [150, 150, 60, 40];
  STRATUM_NOTE = '1 none | 2 dieback 2015-16 | 3 other-year clean | '
               + '4 other-year suspect (2001/2011/2013)';
} else {
  CLASS_VALUES = [1, 2, 3];
  CLASS_POINTS = [150, 150, 100];
  STRATUM_NOTE = '1 none | 2 dieback 2015-16 | 3 other year';
}
print('Stratification:', STRATUM_NOTE);
print('Allocation:', CLASS_POINTS);
 
 
// =====================================================================
// 4 — STRATUM BAND
// =====================================================================

 
var mapClass = classImg.select('class').rename('map_class');
var suspect  = historyImg.select('suspect').unmask(0);
 
var stratum = mapClass.rename('stratum');
if (SPLIT_CLASS3) {
  stratum = stratum.where(mapClass.eq(3).and(suspect.eq(1)), 4);
}
stratum = stratum.updateMask(mapClass.mask()).toByte().rename('stratum');
 
 
// =====================================================================
// 5 — STRATUM WEIGHTS   
// =====================================================================

 
print('STRATUM WEIGHTS — pixel count N_h:', stratum.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(),
  geometry: gulf, scale: 30, maxPixels: 1e13, tileScale: 4
}));
 
print('STRATUM WEIGHTS — area (ha):',
  ee.Image.pixelArea().divide(10000).addBands(stratum).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'stratum'}),
    geometry: gulf, scale: 30, maxPixels: 1e13, tileScale: 4
  }));
 
print('MAP CLASS area (ha) — for the collapsed 3x3 matrix:',
  ee.Image.pixelArea().divide(10000).addBands(mapClass).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'map_class'}),
    geometry: gulf, scale: 30, maxPixels: 1e13, tileScale: 4
  }));
 
 
// =====================================================================
// 6 — SAMPLING IMAGE
// =====================================================================

//
// SENTINEL: -1 means "no detection at this pixel". yod is a year and the
// rest are non-negative for a loss segment, so -1 is unambiguous.
 
var histBands = historyImg
  .select(['yod', 'mag', 'dur', 'preval', 'rate', 'dsnr', 'suspect'],
          ['hist_yod', 'hist_mag', 'hist_dur', 'hist_preval',
           'hist_rate', 'hist_dsnr', 'suspect'])
  .unmask(-1);
 
var dbBands = diebackImg
  .select(['yod', 'mag', 'dur', 'preval', 'rate', 'dsnr'],
          ['db_yod', 'db_mag', 'db_dur', 'db_preval', 'db_rate', 'db_dsnr'])
  .unmask(-1);
 
var sampleImg = stratum
  .addBands(mapClass)
  .addBands(histBands)
  .addBands(dbBands)
  .updateMask(mapClass.mask());
 
 
// =====================================================================
// 7 — DRAW THE SAMPLE
// =====================================================================
var rawSample = sampleImg.stratifiedSample({
  numPoints:   0,              // ignored where classPoints is supplied
  classBand:   'stratum',
  region:      gulf,
  scale:       30,
  seed:        SEED,
  classValues: CLASS_VALUES,
  classPoints: CLASS_POINTS,
  dropNulls:   true,
  tileScale:   16,             
                              
  geometries:  true
});
 
 
// =====================================================================
// 8 — POINT IDS
// =====================================================================
 
var shuffled = rawSample.randomColumn('rnd', SEED).sort('rnd');
var asList   = shuffled.toList(1000);
var n        = shuffled.size();
 
var indexed = ee.FeatureCollection(
  ee.List.sequence(0, n.subtract(1)).map(function (i) {
    i = ee.Number(i);
    var f   = ee.Feature(asList.get(i));
    var xy  = f.geometry().coordinates();
    var lon = ee.Number(xy.get(0));
    var lat = ee.Number(xy.get(1));
    return f
      .set('point_id', i.add(1))
      .set('lon', lon)
      .set('lat', lat)
      // Satellite basemap at the point. Saves pasting coordinates 400
      // times. Google Earth Pro is still where the historical slider is.
      .set('gmaps', ee.String('https://www.google.com/maps/@')
        .cat(lat.format('%.5f')).cat(',')
        .cat(lon.format('%.5f')).cat(',400m/data=!3m1!1e3'));
  })
);
 
 
// =====================================================================
// 9 — ATTACH THE NBR TRAJECTORY
// =====================================================================

 
var sample = nbrStack.reduceRegions({
  collection: indexed,
  reducer:    ee.Reducer.first(),
  scale:      30,
  tileScale:  4
});
 

sample = sample.map(function (f) {
  return f.set({
    ref_class:      '',   
    ref_year:       '',   
    ref_confidence: '',   // 1 low | 2 medium | 3 high
    ref_imagery:    '',   
    ref_notes:      ''
  });
});
 
 
// =====================================================================
// 10 — EXPORTS
// =====================================================================
var CORE = ['point_id', 'lon', 'lat', 'gmaps'];
var REF  = ['ref_class', 'ref_year', 'ref_confidence', 'ref_imagery',
            'ref_notes'];
var MAP  = ['stratum', 'map_class',
            'hist_yod', 'hist_mag', 'hist_dur', 'hist_preval',
            'hist_rate', 'hist_dsnr', 'suspect',
            'db_yod', 'db_mag', 'db_dur', 'db_preval',
            'db_rate', 'db_dsnr'];
 
// FILE A — the blind interpretation sheet
Export.table.toDrive({
  collection:  sample,
  description: 'sample400_blind',
  folder:      'mangrove_resilience',
  fileFormat:  'CSV',
  selectors:   CORE.concat(nbrBandNames).concat(REF)
});
 
// FILE B — the full record
Export.table.toDrive({
  collection:  sample,
  description: 'sample400_full',
  folder:      'mangrove_resilience',
  fileFormat:  'CSV',
  selectors:   CORE.concat(MAP).concat(nbrBandNames)
});
 
// FILE C — asset copy
Export.table.toAsset({
  collection:  sample,
  description: 'sample400_asset',
  assetId:     ASSET_DIR + 'sample400'
});
 
 
// =====================================================================
// 11 — CHECKS BEFORE THE EXPORT TASKS
// =====================================================================
print('Sample size :', rawSample.size());
print('Points per stratum (must match the allocation):',
  rawSample.aggregate_histogram('stratum'));
 

print('db_yod among stratum-2 points :',
  rawSample.filter(ee.Filter.eq('stratum', 2)).aggregate_histogram('db_yod'));
 
// Year spread of the class-3 points
print('hist_yod among class-3 points:',
  rawSample.filter(ee.Filter.eq('map_class', 3))
           .aggregate_histogram('hist_yod'));
 
// Suspect-year share of the class-3 sample
print('Suspect flag among class-3 points:',
  rawSample.filter(ee.Filter.eq('map_class', 3))
           .aggregate_histogram('suspect'));
 
Map.addLayer(classImg, {min: 1, max: 3,
  palette: ['cccccc', 'ff0000', 'ffaa00']}, '3-class map (asset)', true);
 

Map.addLayer(rawSample, {color: '0000ff'}, 'sample 400', false);
 
