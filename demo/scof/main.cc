#include <iostream>
#include <fstream>

#include <SingularityConstrainedOctahedralField.hh>

bool isInFileGood(const std::string& filename)
{
  std::ifstream is(filename.c_str());
  return is.good();
}

bool isOutFileGood(const std::string& filename)
{
  std::ofstream os(filename.c_str());
  return os.good();
}

void printUsage(const std::string& progname)
{
  std::cout << "usage: " << progname << " <inFile> <outFile>" << std::endl
            << std::endl
            << "Reads input tet mesh with properties from <inFile> " << std::endl
            << "and writes the resulting tet mesh with properties to <outFile> in ovm format." << std::endl;
}

int main(int argc, const char* argv[])
{

  if (argc != 3)
  {
    printUsage(argv[0]);
    return 0;
  }

  auto inFilename  = std::string(argv[1]);
  auto outFilename = std::string(argv[2]);

  if (!isInFileGood(inFilename))
  {
    std::cout << "Could not open input File " << inFilename << std::endl;
    return 1;
  }

  if (!isOutFileGood(outFilename))
  {
    std::cout << "Could not open output File " << outFilename << std::endl;
    return 1;
  }

  if (isInFileGood(inFilename) && isOutFileGood(outFilename))
      SCOF::generateSCOF(inFilename, outFilename);

  return 0;
}
