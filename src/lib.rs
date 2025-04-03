use geo::coordinate_position::{CoordPos, CoordinatePosition};
use geo::{coord, Area, BoundingRect, Coord};
use ndarray::{Array2};
use tobj;
use std::collections::{HashSet, HashMap};
use ordered_float::OrderedFloat;
use pyo3::prelude::*;
use numpy::{IntoPyArray, PyArray2};

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

    pub fn new(nrows: usize, ncols: usize, cellsize: f64, nodata: f64, layers: usize) -> Self {
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
}

fn raster_to_python(raster: &Raster) -> PyResult<PyObject> {
    Python::with_gil(|py| {
        // Create a 3D numpy array for the raster layers
        let numpy = py.import("numpy")?;

        // Stack 2D arrays into a 3D numpy array
        let layers: Vec<_> = raster.arrays.iter()
            .map(|array| array.to_owned().into_pyarray(py))
            .collect();

        let stacked_layers = numpy.call_method1("stack", (layers, 0))?;

        // Return the 3D array to Python
        println!("Sending raster layers as a 3D array to Python.");

        Ok(stacked_layers.into_py(py))
    })
}

fn rasterize_faces(raster: &mut Raster, triangles: &Triangles, nodata: f64){
    // Loop over triangulated faces and rasterize them
    let mut height_map: HashMap<(usize, usize), HashSet<OrderedFloat<f64>>> = HashMap::new();

    println!("\nRasterizing faces...");
    for face in 0..triangles.faces.len() {
        let triangle = triangles.get_triangle_geo(face);
        // Get candidate cells from triangle bbox
        let tri_bbox = triangle.bounding_rect();
        let colstart = (tri_bbox.min().x.abs() / raster.cellsize).floor() as usize;
        let colend = (tri_bbox.max().x.abs() / raster.cellsize).ceil() as usize;
        let rowstart = (tri_bbox.min().y.abs() / raster.cellsize).floor() as usize;
        let rowend = (tri_bbox.max().y.abs() / raster.cellsize).ceil() as usize;
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
    }

    for ((i, j), heights) in height_map {
        let mut heights_vec: Vec<_> = heights.into_iter().collect();
        //let mut heights_vec: Vec<_> = heights.into_iter().filter(|&x| x != 0.0).collect();

        // Sort descending to get highest first
        heights_vec.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // Check if more layers are needed
        if heights_vec.len() > raster.layers {
            let new_layers = heights_vec.len();
            raster.arrays.resize(new_layers, Array2::from_elem((raster.nrows, raster.ncols), nodata));
            raster.layers = new_layers;  // Update the layer count
        }

        // Assign highest value to layer 0
        if let Some(highest_val) = heights_vec.get(0) {
            raster.arrays[0][[(raster.nrows - 1 - j), i]] = **highest_val;
        }

        // Assign remaining lower values to subsequent layers (ascending order)
        let mut lower_vals: Vec<_> = heights_vec.iter().skip(1).collect();
        lower_vals.reverse();

        for layer_idx in 1..raster.layers {
            if let Some(val) = lower_vals.get(layer_idx - 1) {
                raster.arrays[layer_idx][[(raster.nrows - 1 - j), i]] = ***val;
            }
        }
    }

}

#[pyfunction]
fn rasterize_from_python(input: String , ncols: usize, nrows: usize, cellsize: f64, origin: [f64; 2], nodataval: f64) -> PyResult<PyObject> {
    let mut raster = Raster {
        ncols,
        nrows,
        cellsize,
        origin,
        nodataval,
        arrays: Vec::new(),
        layers: 0,
    };

    let triangles = load_obj(&input);
    rasterize_faces(&mut raster, &triangles, nodataval);

    raster_to_python(&raster)
}

// Python module entry point
#[pymodule]
fn rusterizer_3d(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rasterize_from_python, m)?)?;
    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use pyo3::prelude::*;
//     use pyo3::types::PyDict;
//
//     #[test]
//     fn test_rasterize_from_python() {
//         pyo3::append_to_inittab!(rusterizer_3d);  // Register the Python module
//         pyo3::prepare_freethreaded_python();      // Initialize the Python interpreter
//
//         Python::with_gil(|py| {
//             let input = "example/flatpoly.obj".to_string();
//             let ncols = 800;
//             let nrows = 800;
//             let cellsize = 0.5;
//             let origin = [0.0, 0.0];
//             let nodataval = -9999.0;
//
//             // Import the 'rusterizer_3d' Python module
//             let rusterizer_3d = py.import("rusterizer_3d")
//                 .expect("Failed to import rusterizer_3d");
//
//             // Get the 'rasterize_from_python' function from the module
//             let rasterize_from_python = rusterizer_3d
//                 .getattr("rasterize_from_python")
//                 .expect("Failed to get rasterize_from_python function");
//
//             // Call the function with the input parameters
//             let result: PyResult<&PyAny> = rasterize_from_python
//                 .call1((input, ncols, nrows, cellsize, origin, nodataval));
//
//             // Handle the result
//             match result {
//                 Ok(py_result) => {
//                     let numpy = py.import("numpy").expect("Failed to import numpy");
//
//                     // Access the shape attribute of the result
//                     let shape = py_result.getattr("shape")
//                         .expect("Failed to get shape attribute");
//
//                     // Print or assert things about the shape (for example)
//                     println!("Shape: {:?}", shape);
//                 },
//                 Err(e) => {
//                     panic!("Error calling rasterize_from_python: {:?}", e);
//                 }
//             }
//         });
//     }
// }
