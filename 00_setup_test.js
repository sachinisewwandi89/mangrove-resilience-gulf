var gulf = ee.Geometry.Rectangle([135.0, -18.0, 141.5, -12.5]);
Map.centerObject(gulf, 6);
Map.addLayer(gulf, {color: 'red'}, 'Study area draft');