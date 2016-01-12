#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/trace.h>

#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

class PoissonSurfaceReconstructrion
{
public:
PoissonSurfaceReconstructrion(std::string fileName);

std::string reconstruct(FT sm_angle = 20.0, FT sm_radius=30.0, FT sm_distance=0.375);

private:
	void setFilesName(std::string fileName);

	std::string outputFileName;
	std::string inputFileName;

};
