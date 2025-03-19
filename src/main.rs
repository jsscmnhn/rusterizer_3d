use clap::Parser;
use geo::coordinate_position::{CoordPos, CoordinatePosition};
use geo::{coord, Area, BoundingRect, Coord};
use indicatif::ProgressBar;
use ndarray::{Array2, Axis};
use std::fs::File;
use std::io::Write;
use std::time::Instant;
use tobj;
use std::collections::{HashSet, HashMap};
use ordered_float::OrderedFloat;

#[cfg(feature = "with_gdal")]
use gdal::DriverManager;

//-- CLI parser
#[derive(Parser)]
#[command(author, version, about, long_about = None, allow_negative_numbers(true))]
struct Cli {
    #[arg(short, long)]
    input: String,
    #[arg(short, long)]
    output: String,
    #[arg()]
    cellsize: f64,
    #[arg()]
    ncols: usize,
    #[arg()]
    nrows: usize,
    #[arg()]
    layers: usize,
    #[arg(short, long, default_value_t = -9999.)]
    nodata: f64,
    #[arg(short, long, default_value_t = 0.)]
    x_transform: f64,
    #[arg(short, long, default_value_t = 0.)]
    y_transform: f64,
}


//-- Types
type Point3 = [f64; 3];
type Point2 = [f64; 2];
type Face = [usize; 3];
type Triangle = [Point3; 3];

//-- Primitives
struct Bbox {
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
}

//-- Basic functions
// Add Point3
fn add_pts3(avar: Point3, bvar: Point3) -> Point3 {
    assert_eq!(avar.len(), bvar.len(), "Trying to add unequal lengths!");
    avar.iter()
        .zip(bvar.iter())
        .map(|(&a, &b)| a + b)
        .collect::<Vec<f64>>()
        .try_into()
        .unwrap()
}

// Linear interpolation within a triangle
pub fn interpolate_linear(triangle: &Triangle, pt: &Coord) -> f64 {
    let a0: f64 = geo::Triangle::new(
        *pt,
        coord! { x: triangle[1][0], y: triangle[1][1] },
        coord! { x: triangle[2][0], y: triangle[2][1] },
    )
    .unsigned_area();
    let a1: f64 = geo::Triangle::new(
        *pt,
        coord! { x: triangle[0][0], y: triangle[0][1] },
        coord! { x: triangle[2][0], y: triangle[2][1] },
    )
    .unsigned_area();
    let a2: f64 = geo::Triangle::new(
        *pt,
        coord! { x: triangle[0][0], y: triangle[0][1] },
        coord! { x: triangle[1][0], y: triangle[1][1] },
    )
    .unsigned_area();

    let mut total = 0.;
    total += triangle[0][2] * a0;
    total += triangle[1][2] * a1;
    total += triangle[2][2] * a2;
    total / (a0 + a1 + a2)
}

//-- Map of OBJ mesh faces
struct Triangles {
    faces: Vec<Face>,
    points: Vec<Point3>,
    bbox: Bbox, // dataset range
}

impl Triangles {
    //-- Constructor using the information from obj import
    pub fn new(firstpt: Point2) -> Self {
        Triangles {
            faces: Vec::new(),
            points: Vec::new(),
            bbox: Bbox {
                xmin: firstpt[0],
                ymin: firstpt[1],
                xmax: firstpt[0],
                ymax: firstpt[1],
            },
        }
    }

    //-- Methods
    // Add face (indices to vertex map self.points)
    pub fn add_face(&mut self, face: Face) {
        self.faces.push(face);
    }

    // Add a vertex
    pub fn add_pt(&mut self, pt: Point3) {
        self.points.push(pt);
        // update bbox
        if pt[0] < self.bbox.xmin {
            self.bbox.xmin = pt[0];
        } else if pt[0] > self.bbox.xmax {
            self.bbox.xmax = pt[0];
        }
        if pt[1] < self.bbox.ymin {
            self.bbox.ymin = pt[1];
        } else if pt[1] > self.bbox.ymax {
            self.bbox.ymax = pt[1];
        }
    }

    // Transform points
    // sets origin of the dataset to [XLL, YLL] and avoids dealing with origin until the output
    pub fn transform_pts(&mut self) {
        let pt_transform = [-self.bbox.xmin, -self.bbox.ymin, 0.];
        for pt in self.points.iter_mut() {
            *pt = *pt // add_pts3(*pt, pt_transform);
        }
    }

    // Return triangle vertices
    // returns 3x3 array [x, y, z] for every face vertex
    pub fn get_triangle(&self, faceidx: usize) -> Triangle {
        assert_eq!(
            self.faces[faceidx].len(),
            3,
            "Triangle structure has more than 3 vertices!"
        );
        [
            self.points[self.faces[faceidx][0]],
            self.points[self.faces[faceidx][1]],
            self.points[self.faces[faceidx][2]],
        ]
    }

    // Return triangle vertices, 2D (x-y) projection in the triangle
    // data struct of package 'geo'
    pub fn get_triangle_geo(&self, faceidx: usize) -> geo::Triangle {
        let pt0 = [
            self.points[self.faces[faceidx][0]][0],
            self.points[self.faces[faceidx][0]][1],
        ];
        let pt1 = [
            self.points[self.faces[faceidx][1]][0],
            self.points[self.faces[faceidx][1]][1],
        ];
        let pt2 = [
            self.points[self.faces[faceidx][2]][0],
            self.points[self.faces[faceidx][2]][1],
        ];
        geo::Triangle::new(
            coord! { x: pt0[0], y: pt0[1] },
            coord! { x: pt1[0], y: pt1[1] },
            coord! { x: pt2[0], y: pt2[1] },
        )
    }

    // Return triangle vertices, height (z-coordinate) only
    /*
    pub fn get_triangle_z(&self, faceidx: usize) -> Vec<f64> {
        let mut triangle: Vec<f64> = Vec::with_capacity(3);
        for ptidx in &self.faces[faceidx] {
            triangle.push(self.points[*ptidx][2]);
        }
        return triangle;
    }
     */
}

//-- Load OBJ into vector of triangles
fn load_obj(filename: &str) -> Triangles {
    println!("Loading file '{}'", filename);
    let load_options = &tobj::LoadOptions {
        triangulate: true,
        ..Default::default()
    };

    let (models, _materials) =
        tobj::load_obj(filename, load_options).expect("Failed to load OBJ file");
    //    println!("Number of models          = {}", models.len());
    let firstpt = &models[0].mesh.positions;
    let mut triangles = Triangles::new([firstpt[0], firstpt[1]]);

    let mut ptstart: usize = 0;
    for (_i, m) in models.iter().enumerate() {
        let mesh = &m.mesh;
        assert_eq!(mesh.indices.len() % 3, 0, "Faces should be triangulated");
        for fidx in 0..mesh.indices.len() / 3 {
            let face_indices: Face = [
                mesh.indices[3 * fidx] as usize + ptstart,
                mesh.indices[3 * fidx + 1] as usize + ptstart,
                mesh.indices[3 * fidx + 2] as usize + ptstart,
            ];
            triangles.add_face(face_indices);
        }
        assert_eq!(
            mesh.positions.len() % 3,
            0,
            "More than three vertices per face!"
        );
        for vtx in 0..mesh.positions.len() / 3 {
            let point = [
                mesh.positions[3 * vtx],
                mesh.positions[3 * vtx + 1],
                mesh.positions[3 * vtx + 2],
            ];
            triangles.add_pt(point);
        }
        ptstart = triangles.points.len();
    }
    // Transform the coordinate system so that origin for calculation
    // is the [XLL, YLL] of the dataset
    triangles.transform_pts();
    return triangles;
}

//-- Raster data structure
struct Raster {
    nrows: usize,
    ncols: usize,
    cellsize: f64,
    origin: Point2,
    nodataval: f64,
    layers: usize,
    arrays: Vec<Array2<f64>>,
}

impl Raster {
    //-- Constructor
    // pub fn new(dataset_range: &Bbox, cellsize: f64, nodata: f64) -> Self {
    //     let nrows = ((dataset_range.ymax - dataset_range.ymin) / cellsize)
    //         .abs()
    //         .ceil() as usize;
    //     let ncols = ((dataset_range.xmax - dataset_range.xmin) / cellsize)
    //         .abs()
    //         .ceil() as usize;
    //     Raster {
    //         nrows,
    //         ncols,
    //         cellsize,
    //         origin: [dataset_range.xmin, dataset_range.ymin],
    //         nodataval: nodata,
    //         array: Array2::from_elem((nrows, ncols), nodata),
    //     }
    // }

    pub fn new(dataset_range: &Bbox, nrows: usize, ncols: usize, cellsize: f64, nodata: f64, layers: usize) -> Self {
        let mut arrays = Vec::new();
        for _ in 0..(layers.max(1)) {
            arrays.push(Array2::from_elem((nrows, ncols), nodata));
        }
        Raster {
            nrows,
            ncols,
            cellsize,
            origin: [0.0, 0.0],
            nodataval: nodata,
            layers,
            arrays,
        }
    }

    //-- Methods
    // Get cell centroid coordinates (x-y) in coord data structure
    // of 'geo' package
    pub fn xy_coord_geo(&self, col: usize, row: usize) -> Coord {
        // Changed: just ignore if a value is out of bounds
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        coord! {
            x: self.cellsize * (0.5 + col as f64),
            y: self.cellsize * (0.5 + row as f64)
        }
    }


    // Write raster to disk in ESRI ASC format
    pub fn write_asc(&self, path: String) -> std::io::Result<()> {
        for (layer_idx, array) in self.arrays.iter().enumerate() {
            let filename = if self.layers > 1 {
                format!(
                    "{}_layer{}.asc",
                    path.trim_end_matches(".asc"),
                    layer_idx
                )
            } else {
                path.clone()
            };
            let mut f = File::create(filename.clone())?;
            let mut s = String::new();
            // Header
            s.push_str(&format!("NCOLS {}\n", self.ncols));
            s.push_str(&format!("NROWS {}\n", self.nrows));
            s.push_str(&format!("XLLCORNER {}\n", self.origin[0]));
            s.push_str(&format!("YLLCORNER {}\n", self.origin[1]));
            s.push_str(&format!("CELLSIZE  {}\n", self.cellsize));
            s.push_str(&format!("NODATA_VALUE {}\n", self.nodataval));
            // Data
            for i in 0..array.dim().0 {
                let col = array
                    .index_axis(Axis(0), i)
                    .iter()
                    .map(|val| format!("{}", val))
                    .collect::<Vec<String>>()
                    .join(" ");
                s.push_str(&format!("{}\n", col));
            }
            // Write to file
            write!(f, "{}", s).unwrap();
            println!("--> Layer {} saved to '{}'", layer_idx, filename);
        }
        Ok(())
    }
}

fn main() {
    let start = Instant::now();
    println!("=== RUSTERIZER ===");

    // Grab input agruments
    let cli = Cli::parse();

    let (input, output, cellsize, ncols, nrows, nodata, layers) = (cli.input, cli.output, cli.cellsize, cli.ncols, cli.nrows, cli.nodata, cli.layers);
    let transform_pt: Point2 = [cli.x_transform, cli.y_transform];

    // Check the output filename
    if !(output.ends_with(".asc") || output.ends_with(".tif")) {
        panic!("Unsupported file format in the output filename! Expected .asc or .tif");
    }

    // Load obj
    let triangles = load_obj(&input);

    // Initialize raster
    let mut raster = Raster::new(&triangles.bbox, nrows, ncols, cellsize, nodata, layers);

    // Print basic info
    println!(
        "Creating a raster of size: [{}, {}]",
        raster.nrows, raster.ncols
    );
    println!(
        "Bbox min: [{}, {}]",
        triangles.bbox.xmin, triangles.bbox.ymin
    );
    println!(
        "Bbox max: [{}, {}]",
        triangles.bbox.xmax, triangles.bbox.ymax
    );
    println!("Number of faces: {:?}", triangles.faces.len());

    // Loop over triangulated faces and rasterize them
    let mut height_map: HashMap<(usize, usize), HashSet<OrderedFloat<f64>>> = HashMap::new();
    let pb = ProgressBar::new(triangles.faces.len() as u64);

    println!("\nRasterizing faces...");
    for face in 0..triangles.faces.len() {
        let triangle = triangles.get_triangle_geo(face);
        // Get candidate cells from triangle bbox
        let tri_bbox = triangle.bounding_rect();
        let colstart = (tri_bbox.min().x.abs() / cellsize).floor() as usize;
        let colend = (tri_bbox.max().x.abs() / cellsize).ceil() as usize;
        let rowstart = (tri_bbox.min().y.abs() / cellsize).floor() as usize;
        let rowend = (tri_bbox.max().y.abs() / cellsize).ceil() as usize;
        //        println!("rowstart - rowend: {} - {}", rowstart, rowend);
        //        println!("colstart - colend: {} - {}", colstart, colend);

        // Check candidate cells
        for i in colstart..colend {
            for j in rowstart..rowend {
                let pt = &raster.xy_coord_geo(i, j);
                let coordpos = triangle.coordinate_position(pt);
                if (coordpos == CoordPos::Inside) || (coordpos == CoordPos::OnBoundary) {
                    // interpolate
                    let height = interpolate_linear(&triangles.get_triangle(face), pt);
                    //                    println!("interpolated height: {} at [{}, {}]", height, i, j);
                    // CHANGED:  assign if the highest value -> save all values
                    height_map.entry((i, j)).or_default().insert(OrderedFloat(height));

                }
            }
        }
        pb.inc(1);
    }
    pb.finish_with_message("done");

    for ((i, j), mut heights) in height_map {
        // descending order
        let mut heights_vec: Vec<_> = heights.into_iter().collect();

        // Sort descending
        heights_vec.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // Assign to raster layers
        for layer_idx in 0..raster.layers.max(1) {
            if let Some(val) = heights_vec.get(layer_idx) {
                raster.arrays[layer_idx][[(raster.nrows - 1 - j), i]] = **val; // Note the double dereference
            }
        }
    }


    // Transform points to the output CRS before writing to disk
    // raster.set_output_origin(transform_pt);

    // Output raster
    println!("\n\nWriting raster to disk...");
    if output.ends_with(".asc") {
        let re = raster.write_asc(output.to_string());
        match re {
            Ok(_x) => println!("--> .asc output saved to '{}'", output),
            Err(_x) => println!("ERROR: path '{}' doesn't exist, abort.", output),
        }
    } else {
        // else it should be .tif as it was checked at the beginning
        // #[cfg(feature = "with_gdal")]
        // {
        //     let re = raster.write_geotiff(output.to_string());
        //     match re {
        //         Ok(_x) => println!("--> .tif output saved to '{}'", output),
        //         Err(_x) => println!("ERROR: path '{}' doesn't exist, abort.", output),
        //     }
        // }
        // #[cfg(not(feature = "with_gdal"))]
        // {
        //     panic!(
        //         "Rusterizer is not compiled with GDAL!\
        //      Use 'cargo build --release --features with_gdal'"
        //     )
        // }
    }

    //    println!("Array: {:?}", raster.array);
    let duration = start.elapsed();
    println!("\nExecution time: {:?}", duration);
    println!("End");
}
