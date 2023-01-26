# Rusterizer

Convert 3D geometries to a 2D raster.

### Data formats
Input: Wavefront OBJ

Output: [Esri ASCII](https://en.wikipedia.org/wiki/Esri_grid)

### Compilation
Clone the repository with submodules
```
git clone --recurse-submodules git@github.com:ipadjen/rusterizer.git
```

To compile, type

```
cargo build --release
```

Access the executable by typing

```
./target/release/rusterizer
```

You can get Cargo [here](https://www.rust-lang.org/tools/install).

### Usage
```rusterizer -i <input> -o <output> <cellsize>```, e.g. ```rusterizer -i example/bk.obj -o bk.asc 0.5```


