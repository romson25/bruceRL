#include "SetPointProcessing.h"

SetPointProcessing::SetPointProcessing(std::string fileName)
{
    setFilesName(fileName);
}
std::string SetPointProcessing::compute(const double outlierRemovedPercentage,
                                        unsigned int outlierNeighborNumber,
                                        const double wlopRetainPercentage,
                                        const double wlopNeighborRadius,
                                        unsigned int smoothNumber,
                                        unsigned int estimateNormalsNumber)
{
    read();
    removeOutlier(outlierRemovedPercentage, outlierNeighborNumber);
    wlopSimplify(wlopRetainPercentage, wlopNeighborRadius);
    smooth(smoothNumber);
    write();

    readPca();
    pcaEstimateNormals(estimateNormalsNumber);
    writePca();

    return outputFileName;
}

void SetPointProcessing::removeOutlier      (const double removedPercentage,
                                             unsigned int neighborNumber)
{
    inPoints.erase(CGAL::remove_outliers
                   (
                       inPoints.begin(), inPoints.end(),
                       CGAL::Identity_property_map<Point>(),
                       neighborNumber, removedPercentage
                    ),
                   inPoints.end());

    std::vector<Point>(inPoints).swap(inPoints);
}
void SetPointProcessing::wlopSimplify       (const double retainPercentage,
                                             const double neighborRadius)
{
    CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>
            (
                inPoints.begin(), inPoints.end(),
                std::back_inserter(outPoints),
                retainPercentage, neighborRadius
            );
}
void SetPointProcessing::smooth             (unsigned int neighborNumber)
{
    CGAL::jet_smooth_point_set<Concurrency_tag>(outPoints.begin(), outPoints.end(), neighborNumber);
}
void SetPointProcessing::pcaEstimateNormals (unsigned int neighborNumber)
{
    CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
                               CGAL::First_of_pair_property_map<PointVectorPair>(),
                               CGAL::Second_of_pair_property_map<PointVectorPair>(),
                               neighborNumber);

    std::list<PointVectorPair>::iterator unoriented_points_begin =
            CGAL::mst_orient_normals(points.begin(), points.end(),
                                     CGAL::First_of_pair_property_map<PointVectorPair>(),
                                     CGAL::Second_of_pair_property_map<PointVectorPair>(),
                                     neighborNumber);

    points.erase(unoriented_points_begin, points.end());
}


void SetPointProcessing::setFilesName   (std::string fileName)
{
    inputFileName       = fileName;
    outputFileName      = fileName;
    temporaryFileName   = fileName;

    const std::string set  = "_set";
    const std::string temp = "_temp";

    int pos = fileName.find_last_of(".");

    outputFileName.insert(pos, set);
    temporaryFileName.insert(pos, temp);
}
void SetPointProcessing::read           ()
{
    std::ifstream stream(inputFileName.c_str());

    if( CGAL::read_xyz_points(stream, std::back_inserter(inPoints)) && !inPoints.empty() )
        std::cout<<"SetPointProcessing::read\tvertices\t"<<inPoints.size()<<std::endl;
    else
        std::cout<<"SetPointProcessing::read\tERROR:\t"<<inputFileName<<std::endl;

    stream.close();
}
void SetPointProcessing::readPca        ()
{
    std::ifstream stream(temporaryFileName.c_str());

    if(CGAL::read_xyz_points(stream, std::back_inserter(points), CGAL::First_of_pair_property_map<PointVectorPair>()) && !points.empty() )
        std::cout<<"SetPointProcessing::readPca\tvertices\t"<<points.size()<<std::endl;
    else
        std::cerr<<"SetPointProcessing::readPca\tERROR:\t"<<temporaryFileName<<std::endl;

    stream.close();
}
void SetPointProcessing::write          ()
{
    std::ofstream stream(temporaryFileName.c_str());

    if(CGAL::write_xyz_points(stream, outPoints.begin(), outPoints.end()))
        std::cout<<"SetPointProcessing::write\tvertices\t"<<outPoints.size()<<std::endl;
    else
        std::cerr<<"SetPointProcessing::write\tERROR:\t"<<temporaryFileName<<std::endl;

    stream.close();
}
std::string SetPointProcessing::writePca()
{
    std::ofstream stream(outputFileName.c_str());

    if(CGAL::write_xyz_points_and_normals
            (
                stream,
                points.begin(),
                points.end(),
                CGAL::First_of_pair_property_map<PointVectorPair>(),
                CGAL::Second_of_pair_property_map<PointVectorPair>()
            )
        )
        std::cout<<"SetPointProcessing::writePca\tvertices\t"<<points.size()<<std::endl;
    else
        std::cerr<<"SetPointProcessing::writePca\tERROR:\t"<<outputFileName<<std::endl;

    stream.close();

return outputFileName;
}
