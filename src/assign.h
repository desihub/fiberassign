// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef ASSIGN_H
#define ASSIGN_H

#include <cstdint>

#include <string>
#include <vector>
#include <map>
#include <memory>

#include <utils.h>
#include <hardware.h>
#include <tiles.h>
#include <targets.h>


namespace fiberassign {

// This class holds the current assignment information and methods for
// refinement.

class Assignment : public std::enable_shared_from_this <Assignment> {

    public :

        typedef std::shared_ptr <Assignment> pshr;

        Assignment(Targets::pshr tgs, TargetsAvailable::pshr tgsavail,
                   LocationsAvailable::pshr locavail);

        void assign_unused(uint8_t tgtype, int32_t max_per_petal = -1,
                           std::string const & pos_type = std::string("POS"),
                           int32_t start_tile = -1, int32_t stop_tile = -1);

        void assign_force(uint8_t tgtype, int32_t required_per_petal = 0,
                          int32_t start_tile = -1, int32_t stop_tile = -1);

        void redistribute_science(int32_t start_tile = -1,
                                  int32_t stop_tile = -1);

        Hardware::pshr hardware() const;

        Targets::pshr targets() const;

        Tiles::pshr tiles() const;

        TargetsAvailable::pshr targets_avail() const;

        LocationsAvailable::pshr locations_avail() const;

        std::vector <int32_t> tiles_assigned() const;

        std::map <int32_t, int64_t> const & tile_location_target(int32_t tile) const;

        std::map < int32_t, std::map <int32_t, int64_t> > loc_target;

        std::map < int64_t, std::map <int32_t, int32_t> > target_loc;

        std::map < int32_t, std::map <int64_t, std::pair <double, double> > >
            tile_target_xy;

    private :

        void parse_tile_range(int32_t start_tile, int32_t stop_tile,
                              int32_t & tstart, int32_t & tstop);

        bool ok_to_assign (Hardware const * hw, int32_t tile, int32_t loc,
            int64_t target, std::map <int64_t,
            std::pair <double, double> > const & target_xy) const;

        int64_t find_best (Hardware const * hw, Targets * tgs,
            int32_t tile, int32_t fiber, uint8_t type,
            std::map <int64_t, std::pair <double, double> > const & target_xy,
            std::vector <int64_t> const & avail) const;

        void assign_tileloc(Hardware const * hw, Targets * tgs, int32_t tile,
            int32_t loc, int64_t target, uint8_t type);

        void unassign_tileloc(Hardware const * hw, Targets * tgs,
            int32_t tile, int32_t loc, uint8_t type);

        void targets_to_project(Targets const * tgs,
            std::map <int32_t, std::vector <int64_t> > const & tgsavail,
            std::vector <int32_t> const & locs,
            std::vector <int64_t> & tgids, std::vector <double> & tgra,
            std::vector <double> & tgdec) const;

        void project_targets(Hardware const * hw,
                             Targets const * tgs,
                             TargetsAvailable const * tgsavail,
                             int32_t tile_id, double tile_ra,
                             double tile_dec,
                             std::map <int64_t, std::pair <double, double> >
                                & target_xy) const;

        void reassign_science_target(int32_t tstart, int32_t tstop,
            int32_t tile, int32_t loc, int64_t target, bool balance_petals,
            std::map <int32_t, std::map <int32_t, bool> > const & done,
            int32_t & best_tile, int32_t & best_loc) const;

        // The number of assigned locations per tile and spectrograph (petal)
        // For each target class.
        std::map <uint8_t, std::map <int32_t, int32_t> > nassign_tile;
        std::map <uint8_t,
            std::map <int32_t, std::map <int32_t, int32_t> > > nassign_petal;

        // shared handle to the hardware configuration.
        Hardware::pshr hw_;

        // shared handle to the target properties.
        Targets::pshr tgs_;

        // shared handle to the tile properties.
        Tiles::pshr tiles_;

        // shared handle to the available targets.
        TargetsAvailable::pshr tgsavail_;

        // shared handle to the available locations for targets.
        LocationsAvailable::pshr locavail_;

};

}
#endif
