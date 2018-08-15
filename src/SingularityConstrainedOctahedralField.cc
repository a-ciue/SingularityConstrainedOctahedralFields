/*
 * SingularityConstrainedOctahedralField.cc
 */

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/Attribs/ColorAttrib.hh>
#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include "SingularityConstrainedOctahedralField.hh"

namespace SCOF {

void generateSCOF(const std::string& inFileName, const std::string& outFileName)
{
	using Vec3d = OpenVolumeMesh::Geometry::Vec3d;
	using TetMesh = OVM::TetrahedralGeometryKernel<Vec3d, OVM::TetrahedralMeshTopologyKernel>;

    OpenVolumeMesh::IO::FileManager fm;
    TetMesh tetmesh;
    fm.readFile(inFileName, tetmesh);

    StopWatch<> time;

    OpenVolumeMesh::StatusAttrib status(tetmesh);
	SCOF::SingularityGraphT<TetMesh> sg(tetmesh);
	MatchingsAndAlignmentExtractionT<TetMesh, TriMesh> maae(tetmesh, sg, status);

    time.start();
    maae.determineMatchings();
    std::cerr << "\nExtract matchings and alignments time:    " <<time.stop()/1000<<" s"<<std::endl;


    fm.writeFile(outFileName, tetmesh);
    time.start();
    std::cerr <<"#####Generate frame field..." << std::endl;
    FrameFieldGeneratorT<TetMesh> ffg(tetmesh);
    ffg.generate_frame_field();
    std::cerr << "\nGenerate frame field time:    " <<time.stop()/1000<<" s"<<std::endl;

    std::string outFileNameNew = outFileName;
    if(replace(outFileNameNew, ".ovm", ".ofs"))
        ffg.save_frames(outFileNameNew);
    else
        std::cerr<<"\nError: failed in saving frames. File name should end with '.ovm'!";
}


bool replace(std::string& str, const std::string& from, const std::string& to)
{
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

} // namespace SCOF
