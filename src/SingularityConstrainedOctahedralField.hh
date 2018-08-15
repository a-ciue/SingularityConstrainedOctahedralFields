/*
 * SingularityConstrainedOctahedralField.hh
 */

#ifndef SINGULARITYCONSTRAINEDOCTAHEDRALFIELD_HH
#define SINGULARITYCONSTRAINEDOCTAHEDRALFIELD_HH

#include <string>

#include "Typedefs.hh"
#include "MatchingsAndAlignmentExtractionT.hh"
#include "FrameFieldGeneratorT.hh"

namespace SCOF
{

void generateSCOF(const std::string& inFileName, const std::string& outFileName);

bool replace(std::string& str, const std::string& from, const std::string& to);


} // namespace SCOF



#endif /* SINGULARITYCONSTRAINEDOCTAHEDRALFIELD_HH */
