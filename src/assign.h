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

// bitmask values returned by check_collisions:
// this positioner is STUCK or BROKEN
const int COLLISION_SELF_STUCK = 1;
// this target would cause a collision with a STUCK or BROKEN positioner
const int COLLISION_WITH_STUCK = 2;
// this target would cause a collision with a petal or GFA edge
const int COLLISION_WITH_EDGES = 4;
// this target would cause a collision with a neighboring not-stuck positioner
const int COLLISION_WITH_NEIGHBOR = 8;

// This class holds the current assignment information and methods for
// refinement.

class Assignment : public std::enable_shared_from_this <Assignment> {

    public :

        typedef std::shared_ptr <Assignment> pshr;

        Assignment(Targets::pshr tgs, TargetsAvailable::pshr tgsavail,
                   LocationsAvailable::pshr locavail,
                   std::map<int32_t, std::map<int32_t, bool> > stuck_sky
                   = std::map<int32_t, std::map<int32_t,bool> >());

        std::map< int32_t, std::map< std::string, int32_t > >
            get_counts(int32_t start_tile = -1, int32_t stop_tile = -1);

        void assign_unused(uint8_t tgtype, int32_t max_per_petal = -1,
                           int32_t max_per_slitblock = -1,
                           std::string const & pos_type = std::string("POS"),
                           int32_t start_tile = -1, int32_t stop_tile = -1,
                           bool use_zero_obsremain = false);

        void assign_force(uint8_t tgtype, int32_t required_per_petal = 0,
                          int32_t required_per_slitblock = 0,
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

        std::map< std::pair< int32_t, int64_t >, int > check_avail_collisions(int32_t tile,
                                                                              bool all_matches=false) const;

        std::map <int32_t, std::map <int32_t, int64_t> > loc_target;

        std::map <int64_t, std::map <int32_t, int32_t> > target_loc;

        std::map <int32_t, std::map <int64_t, std::pair <double, double> > >
            tile_target_xy;

    private :

        typedef std::pair <int32_t, double> location_weight;

        struct location_distance_compare {
            // Define this method here so that it is inline.
            bool operator() (location_weight const & lhs,
                             location_weight const & rhs) const {
                return lhs.second < rhs.second;
            }
        };

        void tile_available(
            int32_t tile_id,
            uint8_t tgtype,
            std::vector <int32_t> const & locs,
            std::map <int32_t, std::vector <target_weight> > & tile_target_avail,
            std::map <int64_t, std::vector <location_weight> > & tile_loc_avail,
            std::vector <target_weight> & tile_target_weights, bool use_zero_obsremain
        ) const;

        int32_t petal_count(
            uint8_t tgtype,
            int32_t tile,
            int32_t petal
        ) const;

        bool petal_count_max(
            uint8_t tgtype,
            int32_t max_per_petal,
            int32_t tile,
            int32_t petal
        ) const;

        int32_t slitblock_count(
            uint8_t tgtype,
            int32_t tile,
            int32_t petal,
            int32_t slitblock
        ) const;

        bool slitblock_count_max(
            uint8_t tgtype,
            int32_t max_per_slitblock,
            int32_t tile,
            int32_t petal,
            int32_t slitblock
        ) const;

        int check_collisions(
            Hardware const * hw,
            int32_t tile,
            int32_t loc,
            int64_t target,
            std::map <int64_t, std::pair <double, double> > const & target_xy,
            bool check_assigned_neighbors = true
        ) const;

        bool ok_to_assign(
            Hardware const * hw,
            int32_t tile,
            int32_t loc,
            int64_t target,
            std::map <int64_t, std::pair <double, double> > const & target_xy
        ) const;

        void assign_tileloc(
            Hardware const * hw,
            Targets * tgs,
            int32_t tile,
            int32_t loc,
            int64_t target,
            uint8_t type
        );

        void unassign_tileloc(
            Hardware const * hw,
            Targets * tgs,
            int32_t tile,
            int32_t loc,
            uint8_t type
        );

        void reassign_science_target(
            int32_t tstart,
            int32_t tstop,
            int32_t tile,
            int32_t loc,
            int64_t target,
            bool force,
            int32_t & new_tile,
            int32_t & new_loc
        ) const;

        // The number of assigned locations per tile and spectrograph (petal)
        // For each target class.
        std::map <uint8_t, std::map <int32_t, int32_t> > nassign_tile;
        // [target_type][tile_id][petal_id] = count
        std::map <uint8_t,
            std::map <int32_t, std::map <int32_t, int32_t> > > nassign_petal;
        // [target_type][tile_id][petal_id][slitblock_id] = count
        std::map <uint8_t,
            std::map <int32_t,
            std::map <int32_t, std::map <int32_t, int32_t> > > > nassign_slitblock;

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
