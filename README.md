# Rusterizer

Convert 3D geometries to a 2D raster.

### Data formats
Input: Wavefront OBJ

Output: [Esri ASCII](https://en.wikipedia.org/wiki/Esri_grid), [GeoTIFF](https://en.wikipedia.org/wiki/GeoTIFF)

### Compilation
Clone the repository and compile
```
cargo build --release
```

Optionally, to export to GeoTIFF format, you have to install GDAL and compile using 'with_gdal' feature
```
cargo build --release --features with_gdal
```


Access the executable by typing
```
./target/release/rusterizer
```

You can get Cargo [here](https://www.rust-lang.org/tools/install).

### Usage
```rusterizer -i <input> -o <output> <cellsize>```, e.g. ```rusterizer -i example/bk.obj -o bk.asc 0.5```
