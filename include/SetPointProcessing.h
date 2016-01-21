#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <boost/property_map/property_map.hpp>

#include <CGAL/remove_outliers.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>
#include <utility>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> PointVectorPair;

    #ifdef CGAL_LINKED_WITH_TBB
    typedef CGAL::Parallel_tag Concurrency_tag;
    #else
    typedef CGAL::Sequential_tag Concurrency_tag;
    #endif

class SetPointProcessing
{
public:
    SetPointProcessing(std::string fileName);

    std::string compute	(const double outlierRemovedPercentage  = 5.0,
                         unsigned int outlierNeighborNumber     = 24,
                         const double wlopRetainPercentage 	= 20.0,
                         const double wlopNeighborRadius	= 0.5,
                         unsigned int smoothNumber		= 24,
                         unsigned int estimateNormalsNumber	= 18,
                         const double sharpnessAngle            = 25,
                         const double edgeSensitivity           = 0,
                         const double neighborRadius            = 0.25,
                         const double inputOutputMultiplier     = 4);

private:
    std::string inputFileName;
    std::string temporaryFileName;
    std::string outputFileName;

    std::vector<Point> inPoints         {};
    std::vector<Point> outPoints        {};
    std::list<PointVectorPair> points	{};

    void removeOutlier      (const double removedPercentage, unsigned int neighborNumber);
    void wlopSimplify       (const double retainPercentage,  const double neighborRadius);
    void smooth             (unsigned int neighborNumber);
    void pcaEstimateNormals (unsigned int neighborNumber);
    void upsampleAndEdge    (const double sharpness_angle,
                             const double edge_sensitivity,
                             const double neighbor_radius,
                             const double input_output_multiplier);

    void setFilesName       (std::string fileName);
    void read();
    void readPca();
    void write();
    std::string writePca();
};
