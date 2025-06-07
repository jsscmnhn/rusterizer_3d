# Rusterizer

Convert 3D geometries to a layered 2d raster, required for the conversion of .obj files in [SOLWEIG_SOLFD](https://github.com/jsscmnhn/SOLWEIG_SOLFD).

Original code by [Ivan Pađen](https://github.com/ipadjen)

### Data formats
Input: Wavefront OBJ

Output: [Esri ASCII](https://en.wikipedia.org/wiki/Esri_grid), [GeoTIFF](https://en.wikipedia.org/wiki/GeoTIFF)

### Compilation
Clone the repository and compile
```
cargo build --release
```

Access the executable by typing
```
./target/release/rusterizer
```

You can get Cargo [here](https://www.rust-lang.org/tools/install).

### Usage
The code is only usefull as python package. A wheel for windows is provided at . To build your own version for MAC or Linux, use (PyO3)[https://github.com/PyO3/maturin]
