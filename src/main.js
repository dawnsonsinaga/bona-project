var START = '2019-01-01';
var END = '2026-10-01';
var CLEAR = 0.65;
var MIN_SHORE = 40;

var aoi = ee.Geometry.Rectangle([98.45, 2.15, 99.20, 2.95]);

var anchors = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([98.5284, 2.8011]), {site: 'Sialahi', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.5235, 2.8672]), {site: 'Paropo', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.5252, 2.8867]), {site: 'Paropo', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.8876, 2.6450]), {site: 'Timur', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.8952, 2.6204]), {site: 'Timur', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.9154, 2.5935]), {site: 'Timur', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.9205, 2.5686]), {site: 'Timur', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([98.9417, 2.5746]), {site: 'Timur', type: 'KJA'}),
  ee.Feature(ee.Geometry.Point([99.0525, 2.3411]), {site: 'K_Balige', type: 'KONTROL'}),
  ee.Feature(ee.Geometry.Point([99.0300, 2.3570]), {site: 'K_BulBul', type: 'KONTROL'}),
  ee.Feature(ee.Geometry.Point([99.0014, 2.3321]), {site: 'K_Meat', type: 'KONTROL'})
]);

var SITES = ['Sialahi', 'Paropo', 'Timur', 'K_Balige', 'K_BulBul', 'K_Meat'];
var TYPE = {
  Sialahi: 'KJA', Paropo: 'KJA', Timur: 'KJA',
  K_Balige: 'KONTROL', K_BulBul: 'KONTROL', K_Meat: 'KONTROL'
};

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

function prepare(img) {
  var sr = img.select(['B3', 'B4', 'B5', 'B6', 'B8', 'B11']).divide(10000);
  var clear = img.select('cs_cdf').gte(CLEAR);
  var water = sr.normalizedDifference(['B3', 'B11']).gt(0);
  var lowNIR = sr.select('B8').lt(0.04);
  var ndci = sr.normalizedDifference(['B5', 'B4']).rename('NDCI');
  var mci = sr.select('B5').subtract(sr.select('B4'))
    .subtract(sr.select('B6').subtract(sr.select('B4')).multiply((705 - 665) / (740 - 665)))
    .multiply(1000)
    .rename('MCI');
  return ndci.addBands(mci)
    .updateMask(clear)
    .updateMask(water)
    .updateMask(lowNIR)
    .copyProperties(img, ['system:time_start']);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(START, END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70))
  .linkCollection(csPlus, ['cs_cdf'])
  .map(prepare);

var occ = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0);
var dShore = occ.lt(80).selfMask()
  .distance(ee.Kernel.euclidean(1500, 'meters'))
  .unmask(1500)
  .rename('dShore');

function siteStats(img, site, label) {
  var own = anchors.filter(ee.Filter.eq('site', site));
  var others = anchors.filter(ee.Filter.and(
    ee.Filter.neq('site', site),
    ee.Filter.eq('type', 'KJA')
  ));
  var d = own.distance(3000);
  var dOther = others.distance(3000).unmask(3000);

  var band = ee.Image(0)
    .where(d.gte(40).and(d.lt(150)), 1)
    .where(d.gte(150).and(d.lt(400)), 2)
    .where(d.gte(400).and(d.lt(800)), 3)
    .where(d.gte(800).and(d.lt(1500)), 4)
    .where(d.gte(1500).and(d.lte(3000)), 5)
    .rename('band');

  var stack = img.select(['NDCI', 'MCI'])
    .addBands(dShore)
    .addBands(band)
    .updateMask(band.gt(0))
    .updateMask(dShore.gte(MIN_SHORE))
    .updateMask(dOther.gte(1500));

  var reducer = ee.Reducer.mean()
    .combine({reducer2: ee.Reducer.count(), sharedInputs: true})
    .repeat(3)
    .group({groupField: 3, groupName: 'band'});

  var res = stack.reduceRegion({
    reducer: reducer,
    geometry: own.geometry().buffer(3000),
    scale: 20,
    maxPixels: 1e10,
    tileScale: 4
  });

  return ee.FeatureCollection(ee.List(res.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var mean = ee.List(g.get('mean'));
    var count = ee.List(g.get('count'));
    return ee.Feature(null, {
      site: site,
      type: TYPE[site],
      periode: label,
      band: g.get('band'),
      NDCI: mean.get(0),
      MCI: mean.get(1),
      jarakTepi_m: mean.get(2),
      n_piksel: count.get(0)
    });
  }));
}

var periods = [{label: 'SEMUA_2019-2026', start: START, end: END}];
for (var y = 2019; y <= 2026; y++) {
  periods.push({
    label: String(y),
    start: y + '-01-01',
    end: y === 2026 ? END : (y + 1) + '-01-01'
  });
}

var collections = [];
periods.forEach(function (p) {
  var comp = s2.filterDate(p.start, p.end).median();
  SITES.forEach(function (s) {
    collections.push(siteStats(comp, s, p.label));
  });
});
var results = ee.FeatureCollection(collections).flatten();

var overall = results.filter(ee.Filter.eq('periode', 'SEMUA_2019-2026'));
print('Halo summary 2019-2026', overall);
print(ui.Chart.feature.groups(overall, 'band', 'NDCI', 'site')
  .setChartType('LineChart')
  .setOptions({
    title: 'NDCI by distance band, KJA vs control (2019-2026)',
    hAxis: {title: 'Distance band'},
    vAxis: {title: 'NDCI'},
    pointSize: 5
  }));

Export.table.toDrive({
  collection: results,
  description: 'BONA_v2_halo',
  fileNamePrefix: 'BONA_v2_halo',
  fileFormat: 'CSV',
  selectors: ['site', 'type', 'periode', 'band', 'NDCI', 'MCI', 'jarakTepi_m', 'n_piksel']
});

Map.setOptions('SATELLITE');
Map.centerObject(anchors, 10);
Map.addLayer(anchors.filter(ee.Filter.eq('type', 'KJA')), {color: 'red'}, 'KJA');
Map.addLayer(anchors.filter(ee.Filter.eq('type', 'KONTROL')), {color: 'cyan'}, 'Control');
Map.addLayer(dShore.lt(MIN_SHORE).selfMask(), {palette: ['orange']}, 'Shore mask', false);
