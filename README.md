# CYGNSS Derived SM Using Machine Learning

Steps followed to retrieve SM from CyGNSS data over a 36 km grid (TxSON) using Artificial Neural Networks are given below. Our CyGNSS-derived SM product is also compared with TxSON in-situ data, SMAP and NASA's L3 SM product.

1. Download CyGNSS data -> get_cygnss_raw_or.sh
2. Extract CyGNSS data over a 36 km region (This can be changed within the file) -> within_36km.py
3. Download NDVI -> MOD13A1 from us.earthdata.nasa (MODIS/TERRA Vegetation Indices 16 day L3 Global 500m SIN Grid V006
4. Convert NDVI .hdf file to .tiff file in ArcGIS
5. Reproject NDVI .tiff file to EPSG:4326 -> reproject_ndvi.ipynb
6. Upscale TxSON (in-situ) data to a 9km grid at an hourly scale using Voronoi Method -> upscaled_TxSON_9km.mat
7. Make the dataset for ML models (Contains code to aggregate txson insitu data to a daily scale) -> make_ml_dataset.ipynb
8. Train ANN on data for all months of the year -> ann_all_months_9km.ipynb
9. Comparison with TxSON, SMAP and NASA's CyGNSS L3 product (This comparison is at a daily 36 km grid resolution) -> l3_smap_comp_all_months.ipynb
