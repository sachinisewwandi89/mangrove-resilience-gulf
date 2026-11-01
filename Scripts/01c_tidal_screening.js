// Tidal_screening  —  PHASE 1 revision
//
// Adds tidal-stage screening by masking inundated pixels per image
// before compositing, using MNDWI as a water proxy.


// ---- 1. Study area and masks 
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);

var gmwUnion = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/union/gmw_v3_mng_union');
var gmwCore = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/core/gmw_v3_mng_core');

var mangroveMask = gmwUnion.gt(0).selfMask().clip(gulf);
var coreMask     = gmwCore.gt(0).selfMask().clip(gulf);


// ---- 2. Cloud mask (unchanged) 
function maskL2(img) {
  var qa = img.select('QA_PIXEL');
  var bad = qa.bitwiseAnd(1 << 0).neq(0)      // fill
      .or(qa.bitwiseAnd(1 << 1).neq(0))       // dilated cloud
      .or(qa.bitwiseAnd(1 << 2).neq(0))       // cirrus
      .or(qa.bitwiseAnd(1 << 3).neq(0))       // cloud
      .or(qa.bitwiseAnd(1 << 4).neq(0));      // cloud shadow
  return img.updateMask(bad.not());
}


// ---- 3. Roy et al. (2016) coefficients
var ROY = {
  green: {slope: 0.8483, itcp: 0.0088},
  nir:   {slope: 0.8462, itcp: 0.0412},
  swir1: {slope: 0.8937, itcp: 0.0254},
  swir2: {slope: 0.9071, itcp: 0.0172}
};

function scaleSR(img) {
  return img.select('SR_B.').multiply(0.0000275).add(-0.2);
}

var WATER_THRESHOLD = 1.0;


// L5 / L7 : green SR_B2, NIR SR_B4, SWIR1 SR_B5, SWIR2 SR_B7
function prepTMETM(img) {
  var sr = scaleSR(img);
  var nbr   = sr.normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR');
  var mndwi = sr.normalizedDifference(['SR_B2', 'SR_B5']).rename('MNDWI');
  var dry   = mndwi.lte(WATER_THRESHOLD);          // keep non-water only
  return ee.Image(nbr.addBands(mndwi).updateMask(dry))
    .copyProperties(img, ['system:time_start']);
}

// L8 / L9 : green SR_B3, NIR SR_B5, SWIR1 SR_B6, SWIR2 SR_B7
function prepOLI(img) {
  var sr = scaleSR(img);
  var green = sr.select('SR_B3')
                .multiply(ROY.green.slope).add(ROY.green.itcp).rename('G');
  var nir   = sr.select('SR_B5')
                .multiply(ROY.nir.slope).add(ROY.nir.itcp).rename('N');
  var swir1 = sr.select('SR_B6')
                .multiply(ROY.swir1.slope).add(ROY.swir1.itcp).rename('S1');
  var swir2 = sr.select('SR_B7')
                .multiply(ROY.swir2.slope).add(ROY.swir2.itcp).rename('S2');
  var stack = green.addBands(nir).addBands(swir1).addBands(swir2);

  var nbr   = stack.normalizedDifference(['N', 'S2']).rename('NBR');
  var mndwi = stack.normalizedDifference(['G', 'S1']).rename('MNDWI');
  var dry   = mndwi.lte(WATER_THRESHOLD);
  return ee.Image(nbr.addBands(mndwi).updateMask(dry))
    .copyProperties(img, ['system:time_start']);
}


// ---- 4. Collections 
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepOLI);
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepOLI);

var allScreened = l5.merge(l7).merge(l8).merge(l9);


// ---- 5. Annual composites, tide-screened 
var years = ee.List.sequence(1987, 2025);

var annualScreened = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    var sub   = allScreened.filterDate(start, end);

    var med  = sub.select('NBR').median().rename('NBR');
    var nObs = sub.select('NBR').count().rename('n_obs').unmask(0);

    return med.addBands(nObs)
      .updateMask(mangroveMask)
      .clip(gulf)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);

print('Screened collection size :', annualScreened.size());


// COMPARISON DIAGNOSTICS

// D3c. Observations per year AFTER screening 

print(ui.Chart.image.series({
  imageCollection: annualScreened.select('n_obs'),
  region: gulf, reducer: ee.Reducer.mean(), scale: 300
}).setOptions({
  title: 'D3c: observations per dry season AFTER tidal screening',
  vAxis: {title: 'n images'}, hAxis: {title: 'Year'}
}));


// ---- D4c. Mean NBR over stable mangrove AFTER screening 
print(ui.Chart.image.series({
  imageCollection: annualScreened.select('NBR').map(function (i) {
    return i.updateMask(coreMask);
  }),
  region: gulf, reducer: ee.Reducer.mean(), scale: 300
}).setOptions({
  title: 'D4c: mean NBR over stable mangrove AFTER tidal screening',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'}
}));


// ---- D5c. Twenty random pixels

var randomPts20 = mangroveMask.rename('m').stratifiedSample({
  numPoints: 20, classBand: 'm', region: gulf,
  scale: 30, geometries: true, seed: 7
});

print(ui.Chart.image.seriesByRegion({
  imageCollection: annualScreened.select('NBR'),
  regions: randomPts20, reducer: ee.Reducer.first(),
  scale: 30, seriesProperty: 'system:index'
}).setOptions({
  title: 'D5c:Twenty pixels AFTER tidal screening',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'},
  lineWidth: 1, pointSize: 0
}));


// ---- D6c. Original pixel, screened 
var pt = ee.Geometry.Point([135.70194, -15.0896]);
print(ui.Chart.image.series({
  imageCollection: annualScreened.select('NBR'),
  region: pt, reducer: ee.Reducer.first(), scale: 30
}).setOptions({
  title: 'D6c: original pixel AFTER tidal screening',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'},
  lineWidth: 2, pointSize: 4
}));

// ---- What is the MNDWI distribution over mangrove? ----
var mndwi2016 = allScreened.filterDate('2016-05-01','2016-10-31')
                  .select('MNDWI').median().updateMask(mangroveMask);
print(ui.Chart.image.histogram({
  image: mndwi2016, region: gulf, scale: 300, maxBuckets: 50
}).setOptions({title: 'MNDWI over mangrove, 2016 dry season'}));

