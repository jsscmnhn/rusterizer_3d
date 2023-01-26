use clap::Parser;
use indicatif::ProgressBar;
use geo::{coord, Coord, BoundingRect, Area};
use geo::coordinate_position::{CoordinatePosition, CoordPos};
use ndarray::{Array2, Axis};
use std::io::Write;
use std::fs::File;
use std::time::Instant;
use tobj_f64;

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
    #[arg(short, long, default_value_t = -9999.)]
    nodata: f64,
    #[arg(short, long, default_value_t = 0.)]
    x_transform: f64,
    #[arg(short, long, default_value_t = 0.)]
    y_transform: f64
}

//-- Types
type Point3   = [f64; 3];
type Point2   = [f64; 2];
type Face     = [usize; 3];
type Triangle = [Point3; 3];

//-- Primitives
struct Bbox {
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64
}

//-- Basic functions
// Add Point3
fn add_pts3(avar: Point3, bvar: Point3) -> Point3 {
    assert_eq!(avar.len(), bvar.len(), "Trying to add unequal lengths!");
    avar.iter().zip(bvar.iter()).map(|(&a, &b)| a + b)
        .collect::<Vec<f64>>()
        .try_into()
        .unwrap()
}

// Linear interpolation within a triangle
pub fn interpolate_linear(triangle: &Triangle, pt: &Coord) -> f64 {
    let a0: f64 = geo::Triangle::new(
        *pt,
        coord!{ x: triangle[1][0], y: triangle[1][1] },
        coord!{ x: triangle[2][0], y: triangle[2][1] }
    ).unsigned_area();
    let a1: f64 = geo::Triangle::new(
        *pt,
        coord!{ x: triangle[0][0], y: triangle[0][1] },
        coord!{ x: triangle[2][0], y: triangle[2][1] }
    ).unsigned_area();
    let a2: f64 = geo::Triangle::new(
        *pt,
        coord!{ x: triangle[0][0], y: triangle[0][1] },
        coord!{ x: triangle[1][0], y: triangle[1][1] }
    ).unsigned_area();

    let mut total = 0.;
    total += triangle[0][2] * a0;
    total += triangle[1][2] * a1;
    total += triangle[2][2] * a2;
    total / (a0 + a1 + a2)
}

//-- Map of OBJ mesh faces
struct Triangles {
    faces : Vec<Face>,
    points: Vec<Point3>,
    bbox  : Bbox         // dataset range
}

impl Triangles {
    //-- Constructor using the information from obj import
    pub fn new(firstpt: Point2) -> Self {
        Triangles {
            faces : Vec::new(),
            points: Vec::new(),
            bbox  : Bbox {
                xmin: firstpt[0],
                ymin: firstpt[1],
                xmax: firstpt[0],
                ymax: firstpt[1]
            }
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
        let pt_transform = [
            -self.bbox.xmin,
            -self.bbox.ymin,
            0.
        ];
        for pt in self.points.iter_mut() {
            *pt = add_pts3(*pt, pt_transform);
        }
    }

    // Return triangle vertices
    // returns 3x3 array [x, y, z] for every face vertex
    pub fn get_triangle(&self, faceidx: usize) -> Triangle {
        assert_eq!(self.faces[faceidx].len(), 3,
                   "Triangle structure has more than 3 vertices!");
        [
            self.points[self.faces[faceidx][0]],
            self.points[self.faces[faceidx][1]],
            self.points[self.faces[faceidx][2]]
        ]
    }

    // Return triangle vertices, 2D (x-y) projection in the triangle
    // data struct of package 'geo'
    pub fn get_triangle_geo(&self, faceidx: usize) -> geo::Triangle {
        let pt0 = [
            self.points[self.faces[faceidx][0]][0],
            self.points[self.faces[faceidx][0]][1]
        ];
        let pt1 = [
            self.points[self.faces[faceidx][1]][0],
            self.points[self.faces[faceidx][1]][1]
        ];
        let pt2 = [
            self.points[self.faces[faceidx][2]][0],
            self.points[self.faces[faceidx][2]][1]
        ];
        geo::Triangle::new(
            coord!{ x: pt0[0], y: pt0[1] },
            coord!{ x: pt1[0], y: pt1[1] },
            coord!{ x: pt2[0], y: pt2[1] }
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
    let load_options = &tobj_f64::LoadOptions {
        triangulate: true,
        ..Default::default()
    };

    let (models, _materials) = tobj_f64::load_obj(filename, load_options)
        .expect("Failed to load OBJ file");
//    println!("Number of models          = {}", models.len());
    let firstpt = &models[0].mesh.positions;
    let mut triangles = Triangles::new([firstpt[0], firstpt[1]]);

    let mut ptstart: usize = 0;
    for (_i, m) in models.iter().enumerate() {
        let mesh = &m.mesh;
        assert_eq!(mesh.indices.len() % 3, 0, "Faces should be triangulated");
        for fidx in 0..mesh.indices.len() / 3 {
            let face_indices: Face = [
                mesh.indices[3 * fidx]     as usize + ptstart,
                mesh.indices[3 * fidx + 1] as usize + ptstart,
                mesh.indices[3 * fidx + 2] as usize + ptstart,
            ];
            triangles.add_face(face_indices);
        }
        assert_eq!(mesh.positions.len() % 3, 0, "More than three vertices per face!");
        for vtx in 0..mesh.positions.len() / 3 {
            let point = [
                mesh.positions[3 * vtx],
                mesh.positions[3 * vtx + 1],
                mesh.positions[3 * vtx + 2]
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
    nrows     : usize,
    ncols     : usize,
    cellsize  : f64,
    origin    : Point2,
    nodataval : f64,
    array     : Array2::<f64>
}

impl Raster {
    //-- Constructor
    pub fn new(dataset_range: &Bbox, cellsize: f64, nodata: f64) -> Self {
        let nrows = ((dataset_range.ymax - dataset_range.ymin) / cellsize).abs().ceil() as usize;
        let ncols = ((dataset_range.xmax - dataset_range.xmin) / cellsize).abs().ceil() as usize;
        Raster {
            nrows,
            ncols,
            cellsize,
            origin    : [dataset_range.xmin, dataset_range.ymin],
            nodataval : nodata,
            array     : Array2::from_elem((nrows, ncols),
                                          nodata)
        }
    }

    //-- Methods
    // Get cell centroid coordinates (x-y) in coord data structure
    // of 'geo' package
    pub fn xy_coord_geo(&self, col: usize, row: usize) -> Coord {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        coord! {
            x: self.cellsize * (0.5 + col as f64),
            y: self.cellsize * (0.5 + row as f64)
        }
    }

    // Set cell value
    pub fn set_val(&mut self, col: usize, row: usize, val: f64) {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        self.array[[(self.nrows - 1 - row), col]] = val;
    }

    // Return cell value
    pub fn at(&self, col: usize, row: usize) -> &f64 {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        &self.array[[(self.nrows - 1 - row), col]]
    }

    // Set the XLL and YLL to user-defined coordinates
    pub fn set_output_origin(&mut self, transform_pt: Point2) {
        self.origin[0] += transform_pt[0];
        self.origin[1] += transform_pt[1];
    }

    // Write raster to disk in ESRI ASC format
    pub fn write_asc(&self, path: String) -> std::io::Result<()> {
        let mut f = File::create(path)?;
        let mut s = String::new();
        // write header
        s.push_str(&format!("NCOLS {}\n",        self.ncols));
        s.push_str(&format!("NROWS {}\n",        self.nrows));
        s.push_str(&format!("XLLCORNER {}\n",    self.origin[0]));
        s.push_str(&format!("YLLCORNER {}\n",    self.origin[1]));
        s.push_str(&format!("CELLSIZE  {}\n",    self.cellsize));
        s.push_str(&format!("NODATA_VALUE {}\n", self.nodataval));
        // write raster data
        for i in 0..self.array.dim().0 {
            let col = self.array.index_axis(Axis(0), i).iter()
                .map(|val| format!("{}", val))
                .collect::<Vec<String>>()
                .join(" ");
            s.push_str(&format!("{}{}", &col, "\n"));
        }
        // output to file
        write!(f, "{}", s).unwrap();
        Ok(())
    }
}

fn main() {
    let start = Instant::now();
    println!("=== RUSTERIZER ===");

    // Grab input agruments
    let cli = Cli::parse();
    let (input, output, cellsize, nodata) = (cli.input, cli.output, cli.cellsize, cli.nodata);
    let transform_pt: Point2 = [cli.x_transform, cli.y_transform];

    // Load obj
    let triangles = load_obj(&input);

    // Initialize raster
    let mut raster = Raster::new(&triangles.bbox, cellsize, nodata);

    // Print basic info
    println!("Creating a raster of size: [{}, {}]", raster.nrows, raster.ncols);
    println!("Bbox min: [{}, {}]", triangles.bbox.xmin, triangles.bbox.ymin);
    println!("Bbox max: [{}, {}]", triangles.bbox.xmax, triangles.bbox.ymax);
    println!("Number of faces: {:?}", triangles.faces.len());

    // Loop over triangulated faces and rasterize them
    let pb = ProgressBar::new(triangles.faces.len() as u64);
    println!("\nRasterizing faces...");
    for face in 0..triangles.faces.len() {
        let triangle = triangles.get_triangle_geo(face);
        // Get candidate cells from triangle bbox
        let tri_bbox = triangle.bounding_rect();
        let colstart = (tri_bbox.min().x.abs() / cellsize).floor() as usize;
        let colend   = (tri_bbox.max().x.abs() / cellsize).ceil()  as usize;
        let rowstart = (tri_bbox.min().y.abs() / cellsize).floor() as usize;
        let rowend   = (tri_bbox.max().y.abs() / cellsize).ceil()  as usize;
//        println!("rowstart - rowend: {} - {}", rowstart, rowend);
//        println!("colstart - colend: {} - {}", colstart, colend);

        // Check candidate cells
        for i in colstart..colend {
            for j in rowstart..rowend {
                let pt = &raster.xy_coord_geo(i, j);
                let coordpos = triangle.coordinate_position(pt);
                if (coordpos == CoordPos::Inside) || (coordpos == CoordPos::OnBoundary) {
                    // interpolate
                    let height = interpolate_linear(
                        &triangles.get_triangle(face),
                        pt
                    );
//                    println!("interpolated height: {} at [{}, {}]", height, i, j);
                    // assign if the highest value
                    if height > *raster.at(i, j) {
                        raster.set_val(i, j, height);
                    }
                }
            }
        }
        pb.inc(1);
    }
    pb.finish_with_message("done");

    // Transform points to the output CRS before writing to disk
    raster.set_output_origin(transform_pt);

    // Output raster
    println!("\n\nWriting raster to disk...");
    let re = raster.write_asc(output.to_string());
    match re {
        Ok(_x) => println!("--> .asc output saved to '{}'", output),
        Err(_x) => println!("ERROR: path '{}' doesn't exist, abort.", output),
    }
//    println!("Array: {:?}", raster.array);
    let duration = start.elapsed();
    println!("\nExecution time: {:?}", duration);
    println!("End");
}