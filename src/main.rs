use tobj;
use ndarray::{Array2, Axis};
use geo::{coord, Coord, BoundingRect, Area};
use geo::coordinate_position::{CoordinatePosition, CoordPos};
use std::fs::File;
use std::io::Write;
use std::time::Instant;
use indicatif::ProgressBar;
use clap::Parser;

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
    nodata: f64
}

//-- Types
type Point3  = [f64; 3];
type Point2  = [f64; 2];
type Face    = [u64; 3];
//type Triangle = [Point3; 3];
type Triangle = Vec<Point3>; //for now

//-- Primitives
struct Bbox {
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64
}

//-- Basic functions
fn add_pts2(avar: Point2, bvar: Point2) -> Point2 {
    assert_eq!(avar.len(), bvar.len(), "Trying to add unequal lengths!");
    avar.iter().zip(bvar.iter()).map(|(&a, &b)| a + b)
        .collect::<Vec<f64>>()
        .try_into()
        .unwrap()
}

/*
 * Using bbox from package 'geo' instead
fn get_bbox(triangle: &Vec<Point2>) -> Bbox {
    assert!(!triangle.is_empty(), "Trying to calculate bbox of an empty triangle!");
    // Start value
    let mut bbox = Bbox {
        xmin: triangle[0][0],
        ymin: triangle[0][1],
        xmax: triangle[0][0],
        ymax: triangle[0][1]
    };
    for pt in triangle {
        if pt[0] < bbox.xmin {
            bbox.xmin = pt[0];
        } else if pt[0] > bbox.xmax {
            bbox.xmax = pt[0];
        }
        if pt[1] < bbox.ymin {
            bbox.ymin = pt[1];
        } else if pt[1] > bbox.ymax {
            bbox.ymax = pt[1];
        }
    }
    return bbox;
}
 */

// Linear interpolation within a triangle
pub fn interpolate_linear(triangle: &Triangle, pt: &Coord) -> f64 {
    let a0: f64 = geo::Triangle::new(
//        coord!{ x: pt[0], y:pt[1] },
        *pt,
        coord!{ x: triangle[1][0], y: triangle[1][1] },
        coord!{ x: triangle[2][0], y: triangle[2][1] }
    ).unsigned_area();
    let a1: f64 = geo::Triangle::new(
//        coord!{ x: pt[0], y:pt[1] },
        *pt,
        coord!{ x: triangle[0][0], y: triangle[0][1] },
        coord!{ x: triangle[2][0], y: triangle[2][1] }
    ).unsigned_area();
    let a2: f64 = geo::Triangle::new(
//        coord!{ x: pt[0], y:pt[1] },
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
    bbox  : Bbox             // dataset range
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

    // Return triangle vertices
    // returns 3x3 array [x, y, z] for every face vertex
    pub fn get_triangle(&self, faceidx: u64) -> Vec<Point3> {
        assert_eq!(self.faces[faceidx as usize].len(), 3 as usize, "Triangle structure has more than 3 vertices!");
        let mut triangle: Vec<Point3> = Vec::with_capacity(3);
        for ptidx in &self.faces[faceidx as usize] {
            let pt = &self.points[*ptidx as usize];
            triangle.push(*pt);
        }
        return triangle;
    }

    // Return triangle vertices, 2D (x-y) projection
    /*
     * Using data structures from package 'geo'
    pub fn get_triangle_xy(&self, faceidx: usize) -> Vec<Point2> {
        let mut triangle: Vec<Point2> = Vec::with_capacity(3);
        for ptidx in &self.faces[faceidx] {
            triangle.push([self.points[*ptidx as usize][0], self.points[*ptidx as usize][1]]);
        }
        return triangle;
    }
     */

    // Return triangle vertices, 2D (x-y) projection in the triangle
    // data struct of package 'geo'
    pub fn get_triangle_geo(&self, faceidx: usize) -> geo::Triangle {
        let pt0 = [
            self.points[self.faces[faceidx][0] as usize][0],
            self.points[self.faces[faceidx][0] as usize][1]
        ];
        let pt1 = [
            self.points[self.faces[faceidx][1] as usize][0],
            self.points[self.faces[faceidx][1] as usize][1]
        ];
        let pt2 = [
            self.points[self.faces[faceidx][2] as usize][0],
            self.points[self.faces[faceidx][2] as usize][1]
        ];
        geo::Triangle::new(
            coord! {x: pt0[0], y: pt0[1]},
            coord! {x: pt1[0], y: pt1[1]},
            coord! {x: pt2[0], y: pt2[1]}
        )
    }

    // Return triangle vertices, height (z-coordinate) only
    /*
    pub fn get_triangle_z(&self, faceidx: u64) -> Vec<f64> {
        let mut triangle: Vec<f64> = Vec::with_capacity(3);
        for ptidx in &self.faces[faceidx as usize] {
            triangle.push(self.points[*ptidx as usize][2]);
        }
        return triangle;
    }
     */
}

//-- Load OBJ into vector of triangles
//todo see what to do about the double precision here
fn load_obj(filename: &str) -> Triangles {
    println!("Loading file '{}'", filename);
    let load_options = &tobj::LoadOptions {
        triangulate: true,
        ..Default::default()
    };

    let (models, _materials) = tobj::load_obj(filename, load_options)
        .expect("Failed to load OBJ file");
//    println!("Number of models          = {}", models.len());
    let firstpt = &models[0].mesh.positions;
    let mut triangles = Triangles::new([firstpt[0] as f64, firstpt[1] as f64]);

    let mut ptstart: usize = 0;
    for (_i, m) in models.iter().enumerate() {
        let mesh = &m.mesh;
        /*
        println!(
            "model[{}].face_count       = {}",
            i,
            mesh.indices.len() / 3
        );
         */
        /*
        let mut next_face = 0;
        for face in 0..mesh.face_arities.len() {
            let end = next_face + mesh.face_arities[face] as usize;

            let face_indices = &mesh.indices[next_face..end];
            println!(" face[{}].indices          = {:?}", face, face_indices);

            next_face = end;
        }
         */
        assert!(mesh.indices.len() % 3 == 0); // faces should be triangulated
        for fidx in 0..mesh.indices.len() / 3 {
//            let face_indices: Vec<u32> = (&mesh.indices[next_face..end]).to_vec();
//            println!(" face[{}].indices          = {:?}", face, face_indices);
            let face_indices: Face = [
                mesh.indices[3 * fidx] as u64     + ptstart as u64,
                mesh.indices[3 * fidx + 1] as u64 + ptstart as u64,
                mesh.indices[3 * fidx + 2] as u64 + ptstart as u64,
            ];
            triangles.add_face(face_indices);
        }
        /*
        println!(
            "model[{}].positions        = {}",
            i,
            mesh.positions.len() / 3
        );
         */
        assert!(mesh.positions.len() % 3 == 0);
        for vtx in 0..mesh.positions.len() / 3 {
            /*
            println!(
                "              position[{}] = ({}, {}, {})",
                vtx,
                mesh.positions[3 * vtx],
                mesh.positions[3 * vtx + 1],
                mesh.positions[3 * vtx + 2]
            );
             */
            let point = [
                mesh.positions[3 * vtx] as f64,
                mesh.positions[3 * vtx + 1] as f64,
                mesh.positions[3 * vtx + 2] as f64
            ];
            triangles.add_pt(point);
        }
        ptstart = triangles.points.len();
//        dbg!(mesh.positions.len(), mesh.indices.len() / 3, triangles.faces.len(), mesh.indices.len(), triangles.points.len(), ptstart);
//        println!("");
    }
//    println!("triangles faces: {:?}", triangles.faces);
//    println!("triangles points: {:?}", triangles.points);
    return triangles;
}

//-- Raster data structure
struct Raster {
    nrows    : u64,
    ncols    : u64,
    cellsize : f64,
    origin   : Point2,
    nodataval: f64,
    array    : Array2::<f64>
}

impl Raster {
    //-- Constructor
    pub fn new(nrows: u64, ncols: u64, cellsize: f64, origin: Point2, nodata: f64) -> Self {
        Raster {
            nrows    : nrows,
            ncols    : ncols,
            cellsize : cellsize,
            origin   : origin,
            nodataval: nodata,
            array    : Array2::from_elem((nrows as usize, ncols as usize), nodata)
        }
    }

    //-- Methods
    // Get cell centroid coordinates (x-y)
    /*
     * Using data structure from package 'geo' instead
    pub fn get_xy_coord(&self, col: u64, row: u64) -> Point2 {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        add_pts2(self.origin,
            [self.cellsize * (0.5 + col as f64), self.cellsize * (0.5 + row as f64)])
    }
     */

    pub fn xy_coord_geo(&self, col: u64, row: u64) -> Coord {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        let pt = add_pts2(self.origin,
            [self.cellsize * (0.5 + col as f64), self.cellsize * (0.5 + row as f64)]);
        coord! { x: pt[0], y: pt[1] }
    }

    // Set cell value
    pub fn set_val(&mut self, col: u64, row: u64, val: f64) {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        self.array[[(self.nrows - 1 - row) as usize, col as usize]] = val;
    }

    // Return cell value
    pub fn at(&self, col: u64, row: u64) -> &f64 {
        assert!(row < self.nrows, "Invalid row index!");
        assert!(col < self.ncols, "Invalid col index!");
        &self.array[[(self.nrows - 1 - row) as usize, col as usize]]
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

    /*
    // Debug hardcoded data
    let input = "testfile.obj";
    let output = "output.asc";
    let cellsize= 0.5;
    let nodata : f64 = -9999.;
     */

    // Load obj
    let triangles = load_obj(&input);

    // Initialize raster
    let ncols = ((triangles.bbox.xmax - triangles.bbox.xmin) / cellsize).abs().ceil() as u64;
    let nrows = ((triangles.bbox.ymax - triangles.bbox.ymin) / cellsize).abs().ceil() as u64;
    let origin = [triangles.bbox.xmin, triangles.bbox.ymin];
    let mut raster = Raster::new(nrows, ncols, cellsize, origin, nodata);

    println!("Creating a raster of size: [{}, {}]", nrows, ncols);
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
        let colstart = ((tri_bbox.min().x - origin[0]).abs() / cellsize).floor() as u64;
        let colend   = ((tri_bbox.max().x - origin[0]).abs() / cellsize).ceil()  as u64;
        let rowstart = ((tri_bbox.min().y - origin[1]).abs() / cellsize).floor() as u64;
        let rowend   = ((tri_bbox.max().y - origin[1]).abs() / cellsize).ceil()  as u64;
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
                        &triangles.get_triangle(face as u64),
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