// 01_annual_composites  —  PHASE 1 
// Annual dry-season NBR series, 1987-2025, Gulf of Carpentaria

// ---- 1. Study area -------------------------------------------------
var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -14.0]);
Map.centerObject(gulf, 7);
 
 
// ---- 2. Mangrove mask ----------------------------------------------
// UNION = every pixel that was mangrove at ANY GMW epoch (1996-2020).
// This is the permissive mask. Using the 2020 layer alone would delete
// the pixels that died in 2015-16 and never recovered.
var gmwUnion = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/union/gmw_v3_mng_union');
 
var gmwCore = ee.Image('projects/earthengine-legacy/assets/projects/'
  + 'sat-io/open-datasets/GMW/core/gmw_v3_mng_core');
 
var mangroveMask = gmwUnion.gt(0).selfMask().clip(gulf);
 
Map.addLayer(mangroveMask, {palette: ['00FF00']}, 'GMW union mask', true);
Map.addLayer(gmwCore.gt(0).selfMask().clip(gulf),
             {palette: ['0000FF']}, 'GMW core (stable) mask', false);
 
 
// ---- 3. Cloud / shadow / fill mask (Landsat Collection 2) ----------
// QA_PIXEL bits: 0 fill, 1 dilated cloud, 2 cirrus, 3 cloud,
//                4 cloud shadow, 5 snow, 6 clear, 7 water
function maskL2(img) {
  var qa = img.select('QA_PIXEL');
  var fill    = qa.bitwiseAnd(1 << 0).neq(0);
  var dilated = qa.bitwiseAnd(1 << 1).neq(0);
  var cirrus  = qa.bitwiseAnd(1 << 2).neq(0);
  var cloud   = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow  = qa.bitwiseAnd(1 << 4).neq(0);
  var bad = fill.or(dilated).or(cirrus).or(cloud).or(shadow);
  return img.updateMask(bad.not());
}

// ---- 4. Roy et al. (2016) harmonisation ----------------------------
// Transforms OLI (L8/L9) reflectance into ETM+ (L7) equivalent space,
// so that L5, L7, L8 and L9 can be treated as one continuous sensor.
// Coefficients are for surface reflectance in 0-1 units, applied AFTER
// the Collection 2 scaling factors. 
//
//   ETM_equivalent = OLI * slope + intercept
//
// Reference: Roy et al. 2016, Remote Sensing of Environment 185:57-70
var ROY = {
  nir:   {slope: 0.8462, itcp: 0.0412},
  swir2: {slope: 0.9071, itcp: 0.0172}
};
 
// Scale Collection 2 DN to surface reflectance (0-1)
function scaleSR(img) {
  return img.select('SR_B.').multiply(0.0000275).add(-0.2);
}
 
// L5 / L7 : NIR = SR_B4, SWIR2 = SR_B7. Already in ETM+ space.
function prepTMETM(img) {
  var sr = scaleSR(img);
  var nbr = sr.normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR');
  return ee.Image(nbr)
    .copyProperties(img, ['system:time_start'])
    .set('sensor', 'TM_ETM');
}
 
// L8 / L9 : NIR = SR_B5, SWIR2 = SR_B7. Harmonise to ETM+ first.
function prepOLI(img) {
  var sr = scaleSR(img);
  var nir = sr.select('SR_B5')
              .multiply(ROY.nir.slope).add(ROY.nir.itcp).rename('NIR');
  var swir2 = sr.select('SR_B7')
              .multiply(ROY.swir2.slope).add(ROY.swir2.itcp).rename('SWIR2');
  var nbr = nir.addBands(swir2)
               .normalizedDifference(['NIR', 'SWIR2']).rename('NBR');
  return ee.Image(nbr)
    .copyProperties(img, ['system:time_start'])
    .set('sensor', 'OLI');
}
 
 
// ---- 5. Build the merged NBR collection ----------------------------
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepTMETM);
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepOLI);
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
           .filterBounds(gulf).map(maskL2).map(prepOLI);
 
var allNBR = l5.merge(l7).merge(l8).merge(l9);
 
// ---- 6. Annual dry-season composites (May - Oct) --------------------
// Each output image has TWO bands:
//   NBR    - median NBR for that dry season
//   n_obs  - how many clear observations went into the median
// n_obs is early-warning band. 
var years = ee.List.sequence(1987, 2025);
 
// Unmasked version
var annualUnmasked = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 5, 1);
    var end   = ee.Date.fromYMD(y, 10, 31);
    var sub   = allNBR.filterDate(start, end);

    var med   = sub.median().rename('NBR');
    var nObs  = sub.count().rename('n_obs').unmask(0);

    return med.addBands(nObs)
      .clip(gulf)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);

// Masked version 
var annual = annualUnmasked.map(function (img) {
  return img.updateMask(mangroveMask);
});
 
// ---- D1. How big is the study population? 

var areaHa = ee.Image.pixelArea().divide(10000).updateMask(mangroveMask);
print('Mangrove area in AOI (ha):', areaHa.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: gulf, scale: 30, maxPixels: 1e13
}));
 
print('Mangrove pixel count:', mangroveMask.reduceRegion({
  reducer: ee.Reducer.count(), geometry: gulf, scale: 30, maxPixels: 1e13
}));
print('Annual collection size (should be 39):', annual.size());

// D2 
var unionBin = gmwUnion.unmask(0).gt(0);
var coreBin  = gmwCore.unmask(0).gt(0);

// Sanity: are union and core actually different layers?
var areaPx = ee.Image.pixelArea().divide(10000);
print('Union area (ha):', areaPx.updateMask(unionBin.selfMask())
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: gulf,
                 scale: 30, maxPixels: 1e13}));
print('Core area (ha):', areaPx.updateMask(coreBin.selfMask())
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: gulf,
                 scale: 30, maxPixels: 1e13}));

// Pixels that were mangrove at some epochs but not all
var changedMask = unionBin.and(coreBin.not()).selfMask();
print('Changed (union minus core, ha):', areaPx.updateMask(changedMask)
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: gulf,
                 scale: 30, maxPixels: 1e13}));

Map.addLayer(changedMask.clip(gulf), {palette: ['ff00ff']},
             'GMW: changed status', true);

// ---- D3. Observations per year -------------------------------------
// Landsat 5 ended Nov 2011, Landsat 8 began Apr 2013 - 2012 is SLC-off L7 only
print(ui.Chart.image.series({
  imageCollection: annual.select('n_obs'),
  region: gulf,
  reducer: ee.Reducer.mean(),
  scale: 300
}).setOptions({
  title: 'D3: mean clear observations per dry season',
  vAxis: {title: 'n images'}, hAxis: {title: 'Year'}
}));

// ---- D4. Sensor-step check -----------------------------------------
// Mean NBR across STABLE mangrove (core mask) for every year. Real
// ecology should be smooth. A jump at 1999 (L5->L7) or 2013 (L7->L8)
// means harmonisation is not working.
var coreMask = gmwCore.gt(0).selfMask().clip(gulf);
print(ui.Chart.image.series({
  imageCollection: annual.select('NBR').map(function (i) {
    return i.updateMask(coreMask);
  }),
  region: gulf,
  reducer: ee.Reducer.mean(),
  scale: 300
}).setOptions({
  title: 'D4: mean NBR over stable mangrove - look for steps at 1999 / 2013',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'}
}));

// ---- D4b. Same chart, UNION mask (all mangrove, including the dead) ----
print(ui.Chart.image.series({
  imageCollection: annual.select('NBR'),
  region: gulf,
  reducer: ee.Reducer.mean(),
  scale: 300
}).setOptions({
  title: 'D4b: mean NBR over ALL mangrove (union)',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'}
}));

// ---- D5. Five RANDOM mangrove pixels ------------------------------

var randomPts = mangroveMask.rename('m').stratifiedSample({
  numPoints: 5,
  classBand: 'm',
  region: gulf,
  scale: 30,
  geometries: true,
  seed: 42
});

print('D5: number of random points:', randomPts.size());
Map.addLayer(randomPts, {color: 'yellow'}, 'D5: random check pixels', true);

print(ui.Chart.image.seriesByRegion({
  imageCollection: annual.select('NBR'),
  regions: randomPts,
  reducer: ee.Reducer.first(),
  scale: 30,
  seriesProperty: 'system:index'
}).setOptions({
  title: 'D5: five random mangrove pixels, 1987-2025',
  vAxis: {title: 'NBR'},
  hAxis: {title: 'Year'},
  lineWidth: 2,
  pointSize: 3
}));
// ---- D5b. Twenty random pixels
// Four pixels is too thin to calibrate on. Twenty shows the envelope:
// how much of the mask is low-NBR, and how common are the single-year
// spikes seen in D5 series 1 (suspected tidal contamination).

var randomPts20 = mangroveMask.rename('m').stratifiedSample({
  numPoints: 20,
  classBand: 'm',
  region: gulf,
  scale: 30,
  geometries: true,
  seed: 7
});

print('D5b: number of points:', randomPts20.size());
Map.addLayer(randomPts20, {color: 'orange'}, 'D5b: 20 random pixels', false);

print(ui.Chart.image.seriesByRegion({
  imageCollection: annual.select('NBR'),
  regions: randomPts20,
  reducer: ee.Reducer.first(),
  scale: 30,
  seriesProperty: 'system:index'
}).setOptions({
  title: 'D5b: twenty random mangrove pixels, 1987-2025',
  vAxis: {title: 'NBR'},
  hAxis: {title: 'Year'},
  lineWidth: 1,
  pointSize: 0
}));

// ---- D6. Original chart pixel, re-checked 
var pt = ee.Geometry.Point([135.70194, -15.0896]);
Map.addLayer(pt, {color: 'red'}, 'D6: original chart pixel', true);

print('D6: pixel in GMW union?', gmwUnion.reduceRegion({
  reducer: ee.Reducer.first(), geometry: pt, scale: 30}));
print('D6: pixel in GMW core (stable)?', gmwCore.reduceRegion({
  reducer: ee.Reducer.first(), geometry: pt, scale: 30}));

print(ui.Chart.image.series({
  imageCollection: annual.select('NBR'),
  region: pt, reducer: ee.Reducer.first(), scale: 30
}).setOptions({
  title: 'D6: original pixel, WITH harmonisation and union mask',
  vAxis: {title: 'NBR'}, hAxis: {title: 'Year'},
  lineWidth: 2, pointSize: 4}));

print(ui.Chart.image.series({
  imageCollection: annual.select('n_obs'),
  region: pt, reducer: ee.Reducer.first(), scale: 30
}).setOptions({
  title: 'D6b: observations per year at this pixel',
  vAxis: {title: 'n images'}, hAxis: {title: 'Year'}}));


// ---- D8. Extent cross-check against Accad et al. (2019) ------------
// Accad maps 156,365 ha of mangrove in the Qld Gulf at 1:5000 from
// 20-30 cm aerial imagery. Catchments inside our box total ~70,000 ha.
// GMW returned 38,136 ha - about half. Shortfall is narrow creek-line
// stands below SAR resolving power.
var qld = ee.Geometry.Rectangle([138.0, -18.0, 141.5, -14.0]);
print('D8: GMW mangrove, Qld portion of box (ha):',
  ee.Image.pixelArea().divide(10000).updateMask(mangroveMask)
    .reduceRegion({reducer: ee.Reducer.sum(), geometry: qld,
                   scale: 30, maxPixels: 1e13}));
