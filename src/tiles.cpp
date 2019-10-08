// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <hardware.h>
#include <tiles.h>

#include <sstream>

namespace fba = fiberassign;


fba::Tiles::Tiles(std::vector <int32_t> ids, std::vector <double> ras,
                  std::vector <double> decs, std::vector <int32_t> obs,
                  std::vector <std::string> timeobs,
                  std::vector <double> thetaobs,
                  std::vector <double> hourangobs) {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    id = ids;
    ra = ras;
    dec = decs;
    obscond = obs;
    obstime = timeobs;
    obstheta = thetaobs;
    obshourang = hourangobs;

    // Construct the mapping of tile ID to position in the given sequence
    for (size_t i = 0; i < id.size(); ++i) {
        order[id[i]] = i;
        if (logger.extra_debug()) {
            logmsg.str("");
            logmsg << "Tiles:  index " << i << " = ID " << id[i];
            logger.debug_tfg(id[i], -1, -1, logmsg.str().c_str());
        }
    }
}
