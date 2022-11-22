# Rusterizer

Convert 3D geometries to a 2D raster.

### Data formats
Input: Wavefront OBJ

Output: ESRI ASC

### Compilation
Download the repository and do

```
cargo build --release
```

Access the executable by typing

```
./target/release/rusterizer
```

You can get Cargo [here](https://www.rust-lang.org/tools/install).

### Usage
```rusterizer -i <input> -o <output> <cellsize>```, e.g. ```rusterizer -i input.obj -o output.obj 0.5```


