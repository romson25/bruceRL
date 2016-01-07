#include "PoissonSurfaceReconstructrion.h"

PoissonSurfaceReconstructrion::PoissonSurfaceReconstructrion(std::string fileName)
{
setFilesName(fileName);
}

void PoissonSurfaceReconstructrion::setFilesName(std::string fileName)
{
    inputFileName = fileName;
    outputFileName = fileName;
    std::string model = "_model.off";
    int modelLen = model.length();

    int pos = fileName.find_last_of("_");
    outputFileName.replace(pos, modelLen, model);
}
std::string PoissonSurfaceReconstructrion::reconstruct()
{
FT sm_angle = 20.0; // Min triangle angle in degrees.
    FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.

    PointList pointsPSR;
    std::ifstream stream(inputFileName.c_str());

    if (CGAL::read_xyz_points_and_normals(
                              stream,
                              std::back_inserter(pointsPSR),
                              CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())
))
    {
std::cout<<"PoissonSurfaceReconstructrion::read\t zakończono pomyślnie"<<std::endl;
    }
else
{
std::cout<<"PoissonSurfaceReconstructrion::read\t ERROR:\t"<<inputFileName<<std::endl;
}

std::cout<<"PoissonSurfaceReconstructrion::0"<<std::endl;

    Poisson_reconstruction_function function(pointsPSR.begin(), pointsPSR.end(),
                                             CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()) );

std::cout<<"PoissonSurfaceReconstructrion::1"<<std::endl;

    if ( ! function.compute_implicit_function() )
        std::cout<<"Problem 2"<<std::endl;

std::cout<<"PoissonSurfaceReconstructrion::2"<<std::endl;

    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(pointsPSR.begin(), pointsPSR.end(),
                                                       6 /* knn = 1 ring */);

std::cout<<"PoissonSurfaceReconstructrion::3"<<std::endl;

    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

std::cout<<"PoissonSurfaceReconstructrion::4"<<std::endl;

    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);

std::cout<<"PoissonSurfaceReconstructrion::5"<<std::endl;

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*average_spacing,  // Max triangle size
                                                        sm_distance*average_spacing); // Approximation error

std::cout<<"PoissonSurfaceReconstructrion::6"<<std::endl;

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

std::cout<<"PoissonSurfaceReconstructrion::7"<<std::endl;

    if(tr.number_of_vertices() == 0)
        std::cout<<"Problem 1"<<std::endl;

std::cout<<"PoissonSurfaceReconstructrion::8"<<std::endl;

    std::ofstream out(outputFileName.c_str());
    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);
    out << output_mesh;

std::cout<<"PoissonSurfaceReconstructrion::9"<<std::endl;

	return outputFileName;
}
