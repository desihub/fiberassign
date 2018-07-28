// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef TILES_H
#define TILES_H

#include <cstdint>
#include <ctime>

#include <string>
#include <vector>
#include <map>
#include <memory>

#include <hardware.h>


namespace fiberassign {


// This class holds vectors of tile information.

class Tiles : public std::enable_shared_from_this <Tiles> {

    public :

        typedef std::shared_ptr <Tiles> pshr;

        Tiles(Hardware::pshr hw, std::vector <int32_t> ids,
            std::vector <double> ras, std::vector <double> decs,
            std::vector <int32_t> obs);

        std::vector <int32_t> id;
        std::vector <double> ra;
        std::vector <double> dec;
        std::vector <int32_t> obscond;

        Hardware::pshr hardware() const;

        std::map <int32_t, size_t> order;

    private :

        Hardware::pshr hw_;

};

}
#endif
