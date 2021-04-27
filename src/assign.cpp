// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <assign.h>

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstdio>

namespace fba = fiberassign;

namespace fbg = fiberassign::geom;


fba::Assignment::Assignment(fba::Targets::pshr tgs,
                            fba::TargetsAvailable::pshr tgsavail,
                            fba::LocationsAvailable::pshr locavail,
                            std::map < int32_t, std::map <int32_t, bool> > stuck_sky) {
    fba::Timer tm;
    tm.start();

    fba::GlobalTimers & gtm = fba::GlobalTimers::get();
    std::ostringstream gtmname;

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    gtmname.str("");
    gtmname << "Assignment ctor: total";
    gtm.start(gtmname.str());

    tgs_ = tgs;
    tgsavail_ = tgsavail;
    locavail_ = locavail;

    // Get the hardware and tile configuration
    tiles_ = tgsavail_->tiles();
    hw_ = tgsavail_->hardware();

    // Initialize assignment counts

    std::vector <uint8_t> tgtypes;
    tgtypes.push_back(TARGET_TYPE_SCIENCE);
    tgtypes.push_back(TARGET_TYPE_STANDARD);
    tgtypes.push_back(TARGET_TYPE_SKY);
    tgtypes.push_back(TARGET_TYPE_SUPPSKY);
    tgtypes.push_back(TARGET_TYPE_SAFE);

    tile_target_xy.clear();

    size_t ntile = tiles_->id.size();
    for (auto const & tp : tgtypes) {
        nassign_tile[tp].clear();
        nassign_petal[tp].clear();
        nassign_slitblock[tp].clear();
        for (size_t t = 0; t < ntile; ++t) {
            int32_t tile_id = tiles_->id[t];
            nassign_tile[tp][tile_id] = 0;
            nassign_petal[tp][tile_id].clear();
            for (int32_t p = 0; p < hw_->npetal; ++p) {
                nassign_petal[tp][tile_id][p] = 0;
            }
            for (int32_t p = 0; p < hw_->npetal; ++p) {
                nassign_slitblock[tp][tile_id][p].clear();
                for (int32_t s = 0; s < hw_->nslitblock; ++s) {
                    nassign_slitblock[tp][tile_id][p][s] = 0;
                }
            }
            tile_target_xy[tile_id].clear();
            if (tp != TARGET_TYPE_SKY)
                continue;
            // for any stuck positioners that land on good sky,
            // increment the counter
            // None on this tile?
            if (stuck_sky.count(tile_id)==0)
                continue;
            for (auto & st : stuck_sky[tile_id]) {
                // st: < loc_id, bool >
                int32_t loc = st.first;
                bool good_sky = st.second;
                if (!good_sky)
                    continue;
                int32_t petal = hw_->loc_petal[loc];
                int32_t slitblock = hw_->loc_slitblock[loc];
                if (slitblock == -1) {
                    // ETC fiber
                    logmsg.str("");
                    logmsg << "tile " << tile_id << " loc " << loc << " petal " << petal << " slitblock " << slitblock << " is type " << hw_->loc_device_type[loc];
                    logger.debug(logmsg.str().c_str());
                    continue;
                }
                nassign_tile .at(tp).at(tile_id)++;
                nassign_petal.at(tp).at(tile_id).at(petal)++;
                nassign_slitblock.at(tp).at(tile_id).at(petal).at(slitblock)++;
                logmsg.str("");
                logmsg << "tile " << tile_id << " loc " << loc
                       << " on petal " << petal << ", slitblock "
                       << slitblock << " is STUCK on a good sky.";
                logger.debug(logmsg.str().c_str());
            }
        }
    }

    loc_target.clear();
    target_loc.clear();

    auto const * ptiles = tiles_.get();
    auto const * phw = hw_.get();
    auto const * ptgs = tgs_.get();
    auto const * ptgsavail = tgsavail_.get();

    gtmname.str("");
    gtmname << "Assignment ctor: project targets";
    gtm.start(gtmname.str());

    #pragma omp parallel for schedule(dynamic) default(none) shared(ntile, ptiles, phw, ptgs, ptgsavail, logmsg, logger)
    for (size_t t = 0; t < ntile; ++t) {
        int32_t tile_id = ptiles->id[t];
        double tile_ra = ptiles->ra[t];
        double tile_dec = ptiles->dec[t];
        double tile_theta = ptiles->obstheta[t];
        std::map <int64_t, std::pair <double, double> > local_xy;
        project_targets(phw, ptgs, ptgsavail, tile_id, tile_ra, tile_dec,
                        tile_theta, local_xy);
        #pragma omp critical
        {
            logmsg.str("");
            logmsg << "projected targets for tile " << tile_id;
            logger.debug(logmsg.str().c_str());
            tile_target_xy.at(tile_id) = local_xy;
        }
    }

    gtm.stop(gtmname.str());

    gtmname.str("");
    gtmname << "Assignment ctor: total";
    gtm.stop(gtmname.str());

    logmsg.str("");
    logmsg << "Assignment constructor project targets";
    tm.stop();
    tm.report(logmsg.str().c_str());
}


std::vector <int32_t> fba::Assignment::tiles_assigned() const {
    std::vector <int32_t> ret;
    for (auto const & it : loc_target) {
        ret.push_back(it.first);
    }
    return ret;
}


std::map <int32_t, int64_t> const & fba::Assignment::tile_location_target(int32_t tile) const {
    return loc_target.at(tile);
}


void fba::Assignment::tile_available(
    int32_t tile_id,
    uint8_t tgtype,
    std::vector <int32_t> const & locs,
    std::map <int32_t, std::vector <target_weight> > & tile_target_avail,
    std::map <int64_t, std::vector <location_weight> > & tile_loc_avail,
    std::vector <target_weight> & tile_target_weights,
    bool use_zero_obsremain
) const {
    // Reset output objects
    tile_target_avail.clear();
    tile_loc_avail.clear();
    tile_target_weights.clear();

    // positioner center locations in curved coordinates
    auto const & loc_pos = hw_->loc_pos_curved_mm;

    auto const & tile_loctg = tgsavail_->data.at(tile_id);
    for (auto const & loc : locs) {
        for (auto const & tgid : tile_loctg.at(loc)) {
            auto const & tg = tgs_->data.at(tgid);
            if ( ! tg.is_type(tgtype)) {
                // This is not the correct target type.
                continue;
            }
            if (
                (tgtype == TARGET_TYPE_SCIENCE)
                && (tg.obsremain <= 0)
                && ! use_zero_obsremain
            ) {
                // Done observing science observations for this target, and we are
                // not considering science targets with zero obs remaining.
                continue;
            }
            // distance from target to positioner
            double dist = fbg::dist(
                loc_pos.at(loc),
                tile_target_xy.at(tile_id).at(tgid)
            );
            double tot_priority = tg.total_priority();
            tile_loc_avail[tgid].push_back(std::make_pair(loc, dist));
            tile_target_avail[loc].push_back(
                std::make_pair(tgid, tot_priority)
            );
            tile_target_weights.push_back(
                std::make_pair(tgid, tot_priority)
            );
        }
    }

    // Sort available locations for each target by distance from shortest to longest
    location_distance_compare loc_dist_comp;
    for (auto & tgloc : tile_loc_avail) {
        std::stable_sort(tgloc.second.begin(), tgloc.second.end(), loc_dist_comp);
    }

    // Sort available targets for each location by total priority
    target_weight_compare tg_comp;
    for (auto & loctg : tile_target_avail) {
        std::stable_sort(loctg.second.begin(), loctg.second.end(), tg_comp);
    }

    return;
}


int32_t fba::Assignment::petal_count(
    uint8_t tgtype,
    int32_t tile,
    int32_t petal
) const {
    int32_t ret = nassign_petal.at(tgtype).at(tile).at(petal);
    // If assigning SUPP_SKY targets, also include the "regular"
    // sky count on this petal and vice-versa.
    if (tgtype == TARGET_TYPE_SUPPSKY) {
        ret += nassign_petal.at(TARGET_TYPE_SKY).at(tile).at(petal);
    }
    if (tgtype == TARGET_TYPE_SKY) {
        ret += nassign_petal.at(TARGET_TYPE_SUPPSKY).at(tile).at(petal);
    }
    return ret;
}


bool fba::Assignment::petal_count_max(
    uint8_t tgtype,
    int32_t max_per_petal,
    int32_t tile,
    int32_t petal
) const {
    // fba::Logger & logger = fba::Logger::get();
    // std::ostringstream logmsg;
    // bool extra_log = logger.extra_debug();

    // Check petal count limits
    int32_t cur_petal = petal_count(tgtype, tile, petal);

    if (cur_petal >= max_per_petal) {
        // Already have enough objects on this petal
        // if (extra_log) {
        //     logmsg.str("");
        //     logmsg << "petal_count_max: tile " << tile
        //         << ", petal " << petal << " has " << cur_petal
        //         << " targets of type " << (int)tgtype
        //         << " (>= " << max_per_petal << ")";
        //     logger.debug_tfg(tile, -1, -1, logmsg.str().c_str());
        // }
        return true;
    } else {
        return false;
    }
}


int32_t fba::Assignment::slitblock_count(
    uint8_t tgtype,
    int32_t tile,
    int32_t petal,
    int32_t slitblock
) const {
    int32_t ret = nassign_slitblock.at(tgtype).at(tile).at(petal).at(slitblock);
    // If assigning SUPP_SKY targets, also include the "regular"
    // sky count on this slitblock and vice-versa.
    if (tgtype == TARGET_TYPE_SUPPSKY) {
        ret += nassign_slitblock.at(TARGET_TYPE_SKY).at(tile).at(petal).at(slitblock);
    }
    if (tgtype == TARGET_TYPE_SKY) {
        ret += nassign_slitblock.at(TARGET_TYPE_SUPPSKY).at(tile).at(petal).at(slitblock);
    }
    return ret;
}


bool fba::Assignment::slitblock_count_max(
    uint8_t tgtype,
    int32_t max_per_slitblock,
    int32_t tile,
    int32_t petal,
    int32_t slitblock
) const {
    // Check slitblock count limits
    int32_t cur_slitblock = slitblock_count(tgtype, tile, petal, slitblock);
    return (cur_slitblock >= max_per_slitblock);
}


void fba::Assignment::assign_unused(uint8_t tgtype, int32_t max_per_petal,
                                    int32_t max_per_slitblock,
                                    std::string const & pos_type,
                                    int32_t start_tile, int32_t stop_tile,
                                    bool use_zero_obsremain) {
    fba::Timer tm;
    tm.start();

    fba::GlobalTimers & gtm = fba::GlobalTimers::get();
    std::ostringstream gtmname;

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    std::string tgstr = fba::target_string(tgtype);

    // Select locations based on positioner type
    auto device_locs = hw_->device_locations(pos_type);

    logmsg.str("");
    logmsg << "assign unused " << tgstr << ":  considering "
        << device_locs.size()
        << " locations of positioner type \"" << pos_type << "\"";
    logger.info(logmsg.str().c_str());

    gtmname.str("");
    gtmname << "unused " << tgstr << ": total";
    gtm.start(gtmname.str());

    if (max_per_petal < 0) {
        // A negative value indicates that there is no limit.
        max_per_petal = 2147483647;
    }
    if (max_per_slitblock < 0) {
        // A negative value indicates that there is no limit.
        max_per_slitblock = 2147483647;
    }

    // Determine our range of tiles
    int32_t tstart;
    int32_t tstop;
    if (start_tile < 0) {
        tstart = 0;
    } else {
        tstart = tiles_->order.at(start_tile);
    }
    if (stop_tile < 0) {
        tstop = tiles_->id.size() - 1;
    } else {
        tstop = tiles_->order.at(stop_tile);
    }

    logmsg.str("");
    logmsg << "assign unused " << tgstr << ":  working on tiles "
        << start_tile << " (index " << tstart << ") to "
        << stop_tile << " (index " << tstop << ")";
    logger.info(logmsg.str().c_str());

    // Per-tile target availability objects (reset for each tile)
    std::map <int32_t, std::vector <target_weight> > tile_target_avail;
    std::map <int64_t, std::vector <location_weight> > tile_loc_avail;
    std::vector <target_weight> tile_target_weights;

    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];
        double tile_ra = tiles_->ra[t];
        double tile_dec = tiles_->dec[t];

        logmsg.str("");
        logmsg << "assign unused " << tgstr << ": working on tile " << tile_id
            << " at RA/DEC = " << tile_ra << " / " << tile_dec;
        logger.debug(logmsg.str().c_str());

        if ((tgsavail_->data.count(tile_id) == 0)
            || (tgsavail_->data.at(tile_id).size() == 0)) {
            // No targets available for the whole tile.
            if (extra_log) {
                logmsg.str("");
                logmsg << "assign unused " << tgstr << ": tile " << tile_id
                    << " at RA/DEC = " << tile_ra << ", " << tile_dec
                    << ": no available targets";
                logger.debug_tfg(tile_id, -1, -1, logmsg.str().c_str());
            }
            continue;
        }

        gtmname.str("");
        gtmname << "unused " << tgstr << ": local tile availability";
        gtm.start(gtmname.str());

        // Compute the locations which are currently unassigned.

        std::vector <int32_t> loc_unassigned;
        for (auto const & loc : device_locs) {
            if ((loc_target[tile_id].count(loc) == 0) ||
                (loc_target[tile_id].at(loc) < 0)) {
                loc_unassigned.push_back(loc);
            }
        }

        logmsg.str("");
        logmsg << "assign unused " << tgstr << ": tile " << tile_id
            << " considering " << loc_unassigned.size() << " unassigned locations";
        logger.debug(logmsg.str().c_str());

        // Available targets for this tile.  We copy this per-tile data so that
        // we can manipulate it and avoid searching over locations that have already
        // been assigned.

        tile_available(
            tile_id,
            tgtype,
            loc_unassigned,
            tile_target_avail,
            tile_loc_avail,
            tile_target_weights,
            use_zero_obsremain
        );

        // Sort targets by total priority from highest to lowest.

        target_weight_compare tg_compare;
        std::stable_sort(
            tile_target_weights.begin(),
            tile_target_weights.end(),
            tg_compare
        );

        logmsg.str("");
        logmsg << "assign unused " << tgstr << ": tile " << tile_id << " has "
            << tile_target_weights.size() << " available targets for these locs";
        logger.debug(logmsg.str().c_str());

        gtm.stop(gtmname.str());

        // Reference to projected target X/Y locations for this tile.
        auto const & target_xy = tile_target_xy.at(tile_id);

        // Locations available to a single target, declared here and reused to avoid
        // repeated memory allocation.
        std::vector <int32_t> loc_avail;

        // Assign targets in priority order to available positioners.

        int32_t nsuccess = 0;

        for (auto const & tgwit : tile_target_weights) {
            // This target ID
            auto const & tgid = tgwit.first;
            // This weight
            auto const & tgweight = tgwit.second;

            // Look at available locations.  These are already sorted from
            // closest to furthest.
            loc_avail.clear();
            for (auto const & locwt : tile_loc_avail.at(tgid)) {
                if ((loc_target[tile_id].count(locwt.first) > 0) &&
                    (loc_target[tile_id].at(locwt.first) >= 0)) {
                    // Already assigned
                    continue;
                }
                loc_avail.push_back(locwt.first);
            }

            // For each available location from closest to furthest...
            for (auto const & loc : loc_avail) {
                // The petal of this location
                int32_t p = hw_->loc_petal.at(loc);
                // Check petal count limits
                if (petal_count_max(tgtype, max_per_petal, tile_id, p)) {
                    continue;
                }

                // The slitblock of this location
                int32_t s = hw_->loc_slitblock.at(loc);
                // Check slitblock count limits
                if (slitblock_count_max(tgtype, max_per_slitblock, tile_id, p, s)) {
                    continue;
                }

                // Can we assign this location to the target?
                if (ok_to_assign(hw_.get(), tile_id, loc, tgid, target_xy)) {
                    // Yes, assign it
                    assign_tileloc(
                        hw_.get(), tgs_.get(), tile_id, loc, tgid, tgtype
                    );
                    nsuccess++;
                } else {
                    // There must be a collision or some other problem.
                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign unused " << tgstr
                            << ": target " << tgid << ", weight = " << tgweight
                            << ": tile " << tile_id << ", loc " << loc
                            << " NOT ok to assign";
                        logger.debug_tfg(tile_id, loc, tgid, logmsg.str().c_str());
                    }
                }
            }
        }
        logmsg.str("");
        logmsg << "assign unused " << tgstr << ": tile " << tile_id
            << " had " << nsuccess << " successful assignments";
        logger.debug(logmsg.str().c_str());
    }

    gtmname.str("");
    gtmname << "unused " << tgstr << ": total";
    gtm.stop(gtmname.str());

    logmsg.str("");
    if ((max_per_petal < 1000000) && (max_per_slitblock < 1000000)) {
        logmsg << "Assign up to " << max_per_petal << " " << tgstr
               << " targets per petal and " << max_per_slitblock
               << " per slitblock to unused locations";
    } else if (max_per_petal < 1000000) {
        logmsg << "Assign up to " << max_per_petal << " " << tgstr
            << " targets per petal to unused locations";
    } else if (max_per_slitblock < 1000000) {
        logmsg << "Assign up to " << max_per_slitblock << " " << tgstr
            << " targets per slitblock to unused locations";
    } else {
        logmsg << "Assign " << tgstr << " targets to unused locations";
    }
    tm.stop();
    tm.report(logmsg.str().c_str());

    return;
}


// In the case where we have fewer or a comparable number of targets as tile / fibers
// to assign, the early tiles will be dominated by the high priority targets and later
// tiles will have the lower priority targets and may be sparsely populated.
//
// When bumping science targets to place sky and standards, this will result in some
// science targets being bumped and left unassigned.  Instead, this function
// pro-actively moves science targets to available future tile / fibers when the
// future petal has fewer science assignments.  This will redistribute the science
// targets more evenly so that placement of standards and sky will require less
// bumping.

void fba::Assignment::redistribute_science(int32_t start_tile,
                                           int32_t stop_tile) {
    fba::Timer tm;
    tm.start();

    fba::GlobalTimers & gtm = fba::GlobalTimers::get();
    std::ostringstream gtmname;

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    gtmname.str("");
    gtmname << "redistribute science: total";
    gtm.start(gtmname.str());

    // Select locations that are science positioners
    auto device_locs = hw_->device_locations("POS");

    // Determine our range of tiles
    int32_t tstart;
    int32_t tstop;
    if (start_tile < 0) {
        tstart = 0;
    } else {
        tstart = tiles_->order.at(start_tile);
    }
    if (stop_tile < 0) {
        tstop = tiles_->id.size() - 1;
    } else {
        tstop = tiles_->order.at(stop_tile);
    }

    logmsg.str("");
    logmsg << "redist:  working on tiles "
        << start_tile << " (index " << tstart << ") to "
        << stop_tile << " (index " << tstop << ")";
    logger.info(logmsg.str().c_str());

    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];
        double tile_ra = tiles_->ra[t];
        double tile_dec = tiles_->dec[t];

        logmsg.str("");
        logmsg << "redist: working on tile " << tile_id
            << " at RA/DEC = " << tile_ra << " / " << tile_dec;
        logger.debug(logmsg.str().c_str());

        if (nassign_tile[TARGET_TYPE_SCIENCE][tile_id] == 0) {
            // Skip tiles that are fully unassigned.
            if (extra_log) {
                logmsg.str("");
                logmsg << "redist: tile " << tile_id
                    << " at RA/DEC = " << tile_ra << ", " << tile_dec
                    << ": no available targets";
                logger.debug_tfg(tile_id, -1, -1, logmsg.str().c_str());
            }
            continue;
        }

        // Compute the weighting of the assigned science targets.

        std::vector <target_weight> science_targets;

        for (auto const & loc : device_locs) {
            if ((loc_target[tile_id].count(loc) > 0) &&
                (loc_target[tile_id].at(loc) >= 0)) {
                // We have something assigned here...
                auto tgid = loc_target[tile_id].at(loc);
                auto const & tg = tgs_->data.at(tgid);
                if (tg.is_science() && (! tg.is_standard())) {
                    // This is a science target and NOT a standard (we don't
                    // try reassign dual targets)
                    science_targets.push_back(
                        std::make_pair(tgid, tg.total_priority())
                    );
                }
            }
        }

        // Sort the currently assigned science targets by inverse priority order.
        target_weight_rcompare tg_inverse_comp;
        std::stable_sort(
            science_targets.begin(),
            science_targets.end(),
            tg_inverse_comp
        );

        for (auto const & tgwit : science_targets) {
            // This current science target ID
            auto const & tgid = tgwit.first;
            // This weight
            // auto const & tgweight = tgwit.second;

            // Find the location where this is currently assigned
            int32_t tgloc = target_loc.at(tgid).at(tile_id);

            // Try to assign this science target to a future tile
            int32_t new_tile;
            int32_t new_loc;
            reassign_science_target(
                t + 1, tstop, tile_id, tgloc, tgid, false, new_tile, new_loc
            );
            if (new_tile != tile_id) {
                // Some future tile has a better location.  Reassign.
                unassign_tileloc(
                    hw_.get(),
                    tgs_.get(),
                    tile_id,
                    tgloc,
                    TARGET_TYPE_SCIENCE
                );
                assign_tileloc(
                    hw_.get(),
                    tgs_.get(),
                    new_tile,
                    new_loc,
                    tgid,
                    TARGET_TYPE_SCIENCE
                );
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "redist: tile " << tile_id
                        << " loc " << tgloc
                        << " moved science " << tgid
                        << " to tile " << new_tile
                        << ", loc " << new_loc;
                    logger.debug_tfg(
                        tile_id, tgloc, tgid, logmsg.str().c_str()
                    );
                    logger.debug_tfg(
                        new_tile, new_loc, tgid, logmsg.str().c_str()
                    );
                }
            }
        }
    }

    gtmname.str("");
    gtmname << "redistribute science: total";
    gtm.stop(gtmname.str());

    tm.stop();
    tm.report("Redistribute science targets");

    return;
}


void fba::Assignment::assign_force(uint8_t tgtype, int32_t required_per_petal,
                                   int32_t required_per_slitblock,
                                   int32_t start_tile, int32_t stop_tile) {
    fba::Timer tm;
    tm.start();

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    fba::GlobalTimers & gtm = fba::GlobalTimers::get();
    std::ostringstream gtmname;

    std::string tgstr = fba::target_string(tgtype);

    gtmname.str("");
    gtmname << "force " << tgstr << ": total";
    gtm.start(gtmname.str());

    // Select locations that are science positioners
    auto device_locs = hw_->device_locations("POS");

    // Determine our range of tiles
    int32_t tstart;
    int32_t tstop;
    if (start_tile < 0) {
        tstart = 0;
    } else {
        tstart = tiles_->order.at(start_tile);
    }
    if (stop_tile < 0) {
        tstop = tiles_->id.size() - 1;
    } else {
        tstop = tiles_->order.at(stop_tile);
    }

    logmsg.str("");
    logmsg << "assign force " << tgstr << ":  working on tiles "
        << start_tile << " (index " << tstart << ") to "
        << stop_tile << " (index " << tstop << ")";
    logger.info(logmsg.str().c_str());

    // Per-tile target availability objects (reset for each tile)
    std::map <int32_t, std::vector <target_weight> > tile_target_avail;
    std::map <int64_t, std::vector <location_weight> > tile_loc_avail;
    std::vector <target_weight> tile_target_weights;

    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];
        double tile_ra = tiles_->ra[t];
        double tile_dec = tiles_->dec[t];

        logmsg.str("");
        logmsg << "assign force " << tgstr << ": working on tile " << tile_id
            << " at RA/DEC = " << tile_ra << " / " << tile_dec;
        logger.debug(logmsg.str().c_str());

        if (nassign_tile[TARGET_TYPE_SCIENCE][tile_id] == 0) {
            // Skip tiles that are fully unassigned.
            if (extra_log) {
                logmsg.str("");
                logmsg << "assign force " << tgstr << ": tile " << tile_id
                    << " at RA/DEC = " << tile_ra << ", " << tile_dec
                    << ": no available targets";
                logger.debug_tfg(tile_id, -1, -1, logmsg.str().c_str());
            }
            continue;
        }

        gtmname.str("");
        gtmname << "force " << tgstr << ": local tile compute science assignment";
        gtm.start(gtmname.str());

        // Compute the locations which are currently assigned to science targets.
        // Also compute the inverse-priority weighting of these assigned targets.

        std::vector <target_weight> science_targets;
        std::vector <int32_t> loc_science;

        for (auto const & loc : device_locs) {
            if ((loc_target[tile_id].count(loc) > 0) &&
                (loc_target[tile_id].at(loc) >= 0)) {
                // We have something assigned here...
                auto tgid = loc_target[tile_id].at(loc);
                auto const & tg = tgs_->data.at(tgid);
                if (tg.is_science() && (! tg.is_standard())) {
                    // This is a science target and NOT a standard (we don't
                    // try to bump dual targets)
                    loc_science.push_back(loc);
                    science_targets.push_back(
                        std::make_pair(tgid, tg.total_priority())
                    );
                }
            }
        }

        // Sort the currently assigned science targets by inverse priority order
        target_weight_rcompare tg_inverse_comp;
        std::stable_sort(
            science_targets.begin(),
            science_targets.end(),
            tg_inverse_comp
        );

        gtm.stop(gtmname.str());

        gtmname.str("");
        gtmname << "force " << tgstr << ": local tile availability";
        gtm.start(gtmname.str());

        // Available targets for this tile.

        tile_available(
            tile_id,
            tgtype,
            loc_science,
            tile_target_avail,
            tile_loc_avail,
            tile_target_weights,
            false
        );

        gtm.stop(gtmname.str());

        logmsg.str("");
        logmsg << "assign force " << tgstr << ": tile " << tile_id
            << " has " << science_targets.size()
            << " currently assigned science locations";
        logger.debug(logmsg.str().c_str());

        logmsg.str("");
        logmsg << "assign force " << tgstr << ": tile " << tile_id
            << " these locations have " << tile_target_weights.size()
            << " available targets";
        logger.debug(logmsg.str().c_str());

        // Reference to projected target X/Y locations for this tile.
        auto const & target_xy = tile_target_xy.at(tile_id);

        // Targets available to a single location, declared here to avoid
        // repeated memory allocation.
        std::vector <int64_t> tg_avail;

        // The unique petals used for this tile
        std::set <int32_t> unique_petal;
        // The unique (petal,slitblocks) used for this tile
        std::set <std::pair<int32_t, int32_t> > unique_slitblock;

        for (auto const & tgwit : science_targets) {
            // This current science target ID
            auto const & tgid = tgwit.first;
            // This weight
            auto const & tgweight = tgwit.second;

            // Find the location where this is currently assigned
            int32_t tgloc = target_loc.at(tgid).at(tile_id);

            // The petal of this location
            int32_t p = hw_->loc_petal.at(tgloc);
            unique_petal.insert(p);

            // Check petal count limits
            if (petal_count_max(tgtype, required_per_petal, tile_id, p)) {
                continue;
            }

            // The slitblock of this location
            int32_t s = hw_->loc_slitblock.at(tgloc);
            unique_slitblock.insert(std::make_pair(p, s));
            // Check slitblock count limits
            if (slitblock_count_max(tgtype, required_per_slitblock, tile_id, p, s)) {
                continue;
            }

            if (tile_target_avail.count(tgloc) == 0) {
                // There are no available targets of our requested type at this loc.
                continue;
            }

            // Look at available targets for this location.  These are already sorted
            // by total priority.
            tg_avail.clear();
            for (auto const & tgwt : tile_target_avail.at(tgloc)) {
                tg_avail.push_back(tgwt.first);
            }

            // For each available target at this current location...
            for (auto const & avtg : tg_avail) {
                // Can we assign this target?
                if (ok_to_assign(hw_.get(), tile_id, tgloc, avtg, target_xy)) {
                    // Yes, try to assign the bumped science target to a future tile
                    int32_t new_tile;
                    int32_t new_loc;
                    reassign_science_target(
                        t + 1, tstop, tile_id, tgloc, tgid, true, new_tile, new_loc
                    );
                    // Now unassign science target.
                    unassign_tileloc(
                        hw_.get(),
                        tgs_.get(),
                        tile_id,
                        tgloc,
                        TARGET_TYPE_SCIENCE
                    );
                    // Assign the new target
                    assign_tileloc(
                        hw_.get(), tgs_.get(), tile_id, tgloc, avtg, tgtype
                    );
                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign force " << tgstr
                            << ": tile " << tile_id
                            << " loc " << tgloc
                            << " petal " << p
                            << " slitblock " << s
                            << " bumped science " << tgid
                            << " with weight " << tgweight
                            << ", replaced with " << avtg;
                        logger.debug_tfg(tile_id, tgloc, tgid, logmsg.str().c_str());
                        logger.debug_tfg(tile_id, tgloc, avtg, logmsg.str().c_str());
                    }
                    // If we were able, reassign the science target
                    if (new_tile >= 0) {
                        // We were able to find a spot
                        assign_tileloc(
                            hw_.get(),
                            tgs_.get(),
                            new_tile,
                            new_loc,
                            tgid,
                            TARGET_TYPE_SCIENCE
                        );
                        if (extra_log) {
                            logmsg.str("");
                            logmsg << "assign force " << tgstr
                                << ": tile " << tile_id
                                << " loc " << tgloc
                                << " petal " << p
                                << " slitblock " << s
                                << " reassign bumped science " << tgid
                                << " to tile " << new_tile
                                << ", loc " << new_loc;
                            logger.debug_tfg(
                                tile_id, tgloc, tgid, logmsg.str().c_str()
                            );
                            logger.debug_tfg(
                                new_tile, new_loc, tgid, logmsg.str().c_str()
                            );
                        }
                    } else {
                        if (extra_log) {
                            logmsg.str("");
                            logmsg << "assign force " << tgstr
                                << ": tile " << tile_id
                                << " loc " << tgloc
                                << " petal " << p
                                << " slitblock " << s
                                << " bumped science " << tgid
                                << " cannot be reassigned.";
                            logger.debug_tfg(
                                tile_id, tgloc, tgid, logmsg.str().c_str()
                            );
                        }
                    }
                    break;
                } else {
                    // There must be a collision or some other problem.
                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign force " << tgstr
                            << ": tile " << tile_id
                            << " loc " << tgloc
                            << " petal " << p
                            << " slitblock " << s
                            << " cannot bump science " << tgid
                            << " (weight " << tgweight << ")"
                            << " with " << avtg << ": not ok to assign";
                        logger.debug_tfg(tile_id, tgloc, tgid, logmsg.str().c_str());
                        logger.debug_tfg(tile_id, tgloc, avtg, logmsg.str().c_str());
                    }
                }
            }
        }

        // Check final petal assignment counts
        for (auto const & p : unique_petal) {
            int32_t cur_petal = petal_count(tgtype, tile_id, p);
            if (cur_petal < required_per_petal) {
                logmsg.str("");
                logmsg << "assign force " << tgstr << ": tile " << tile_id
                    << ", petal " << p << " could only assign " << cur_petal
                    << " (require " << required_per_petal
                    << ").  Insufficient number of objects or too many collisions";
                logger.warning(logmsg.str().c_str());
            }
        }
        // Check final slitblock assignment counts
        for (auto const & u : unique_slitblock) {
            int32_t p = u.first;
            int32_t s = u.second;
            int32_t cur_slitblock = slitblock_count(tgtype, tile_id, p, s);
            if (cur_slitblock < required_per_slitblock) {
                logmsg.str("");
                logmsg << "assign force " << tgstr << ": tile " << tile_id
                       << ", petal " << p
                       << ", slitblock " << s << " could only assign " << cur_slitblock
                       << " (require " << required_per_slitblock
                       << ").  Insufficient number of objects or too many collisions";
                logger.warning(logmsg.str().c_str());
            }
        }
    }

    gtmname.str("");
    gtmname << "force " << tgstr << ": total";
    gtm.stop(gtmname.str());

    logmsg.str("");
    if (required_per_petal && required_per_slitblock) {
        logmsg << "Force assignment of " << required_per_petal << " "
               << tgstr << " targets per petal and " << required_per_slitblock
               << " per slitblock";
    } else if (required_per_petal) {
        logmsg << "Force assignment of " << required_per_petal << " "
               << tgstr << " targets per petal";
    } else if (required_per_slitblock) {
        logmsg << "Force assignment of " << required_per_slitblock << " "
               << tgstr << " targets per slitblock";
    }

    tm.stop();
    tm.report(logmsg.str().c_str());

    return;
}


void fba::Assignment::reassign_science_target(int32_t tstart, int32_t tstop,
    int32_t tile, int32_t loc, int64_t target, bool force,
    int32_t & new_tile, int32_t & new_loc) const {

    // The "force" option controls whether we want to reassign the science target to
    // a new location even if that location is not as "good" as current assignment.
    // This can happen if the target is being bumped.  If force==true and no new
    // assignment is possible, negative values are returned for new_tile and new_loc.

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    if (extra_log) {
        logmsg.str("");
        logmsg << "reassign: tile " << tile << ", location "
            << loc << ", target " << target << " considering tile indices "
            << tstart << " to " << tstop;
        logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
    }

    // Get the number of unused locations on this current petal.
    int32_t petal = hw_->loc_petal.at(loc);
    int32_t passign = petal_count(TARGET_TYPE_SCIENCE, tile, petal);

    // Vector of available tile / loc pairs which have a loc that is a
    // science positioner.  The available tile / location pairs are already
    // sorted by tile order.
    std::vector < std::pair <int32_t, int32_t> > avail;

    auto const & locavailtg = locavail_->data.at(target);
    std::string pos_str("POS");

    for (auto const & av : locavailtg) {
        int32_t av_tile = av.first;
        int32_t av_tile_indx = tiles_->order.at(av_tile);
        int32_t av_loc = av.second;
        if (nassign_tile.at(TARGET_TYPE_SCIENCE).at(av_tile) == 0) {
            // This available tile / loc is on a tile with
            // nothing assigned.  Skip it.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", location "
                    << loc << ", target " << target
                    << " available tile " << av_tile
                    << " has nothing assigned- skipping";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            continue;
        }
        if (loc_target.at(av_tile).count(av_loc) > 0) {
            // This available tile / loc is already assigned.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " avail tile/loc " << av_tile << "," << av_loc
                    << " already assigned";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            continue;
        }
        if (target_loc.at(target).count(av_tile) > 0) {
            // We have already assigned a location on this tile to this
            // target.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " already assigned on available tile " << av_tile;
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            continue;
        }
        if (av_tile_indx <= tstart) {
            // This available tile came before our current tile.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " available tile " << av_tile
                    << " at index " << av_tile_indx
                    << " is prior to tile start index (" << tstart << ")";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            continue;
        }
        if (pos_str.compare(hw_->loc_device_type.at(av_loc)) == 0) {
            // Available location is not a science positioner
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " available tile " << av_tile
                    << " at index " << av_tile_indx
                    << " is not a science positioner (POS)";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            continue;
        }
        avail.push_back(av);
    }

    new_tile = -1;
    new_loc = -1;
    int32_t best_passign = 500;

    for (auto const & av : avail) {
        // The available tile loc pair
        int32_t av_tile = av.first;
        int32_t av_loc = av.second;
        int32_t av_petal = hw_->loc_petal.at(av_loc);
        // int32_t av_tile_indx = tiles_->order.at(av_tile);

        // Projected target locations on the available tile.
        auto const & av_target_xy = tile_target_xy.at(av_tile);

        if ( ! ok_to_assign(hw_.get(), av_tile, av_loc, target,
                            av_target_xy)) {
            // There must be a collision or some other problem.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " avail tile/loc " << av_tile << "," << av_loc
                    << " not OK to assign";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            continue;
        }

        // At this point we know we have a new tile / loc where we can place
        // this target.  Get the number of assigned locations on the petal of this
        // available tile/loc.
        int32_t av_passign =
            nassign_petal.at(TARGET_TYPE_SCIENCE).at(av_tile).at(av_petal);

        if ((av_passign < hw_->nfiber_petal) && (av_passign < best_passign)) {
            // There are some unassigned locs on this available petal,
            // and the number of unassigned is greater than our current best
            // tile/loc.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " avail tile/loc " << av_tile << "," << av_loc
                    << " new best alternate location for petal counts ("
                    << av_passign << " < " << best_passign << ")";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            new_tile = av_tile;
            new_loc = av_loc;
            best_passign = av_passign;
        } else {
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " avail tile/loc " << av_tile << "," << av_loc
                    << " skipping alternate loc with more petal counts ("
                    << av_passign << " >= " << best_passign << ")";
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
        }
    }

    // If not forcing assignment and we have nothing better, return the original.
    if ((! force) && (passign <= best_passign)) {
        new_tile = tile;
        new_loc = loc;
    }

    return;
}


bool fba::Assignment::ok_to_assign (fba::Hardware const * hw, int32_t tile,
    int32_t loc, int64_t target,
    std::map <int64_t, std::pair <double, double> > const & target_xy
    ) const {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    // Is the location stuck or broken?
    if (
        (hw->state.at(loc) & FIBER_STATE_STUCK) ||
        (hw->state.at(loc) & FIBER_STATE_BROKEN)
    ) {
        if (extra_log) {
            logmsg.str("");
            logmsg << "ok_to_assign: tile " << tile << ", loc "
                << loc << " not OK";
            logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
        }
        return false;
    }

    // NOTE:  When building the TargetsAvailable instance, we already check that the
    // positioner can physically reach every available target.  No need to check that
    // here.

    // Const reference to the assignment for this tile.
    auto const & ftile = loc_target.at(tile);

    std::vector <int32_t> nbs;
    std::vector <int64_t> nbtarget;

    auto const & neighbors = hw->neighbors.at(loc);

    // Check neighboring assignment.  If one of the neighboring positioners is stuck or
    // broken, we keep it for consideration below when checking for collisions.
    for (auto const & nb : neighbors) {
        if (
            (hw->state.at(nb) & FIBER_STATE_STUCK) ||
            (hw->state.at(nb) & FIBER_STATE_BROKEN)
        ) {
            // Include this neighbor in the list to check
            nbs.push_back(nb);
            nbtarget.push_back(-1);
        } else if (ftile.count(nb) > 0) {
            // This neighbor has some assignment.
            int64_t nbtg = ftile.at(nb);
            if (nbtg == target) {
                // Target already assigned to a neighbor.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "ok_to_assign: tile " << tile << ", loc "
                        << loc << ", target " << target
                        << " already assigned to neighbor loc " << nb;
                    logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
                }
                return false;
            }
            nbs.push_back(nb);
            nbtarget.push_back(nbtg);
        }
    }

    // Would assigning this target produce a collision?

    size_t nnb = nbs.size();

    bool collide = false;

    fbg::dpair tpos = target_xy.at(target);

    // On average, the number of neighbors is 2-3.  Threading overhead seems
    // to negate the benefit here.
    // #pragma omp parallel for reduction(||:collide) schedule(dynamic) default(none) shared(hw, nbs, nbtarget, nnb, tcenter, tpos, target_xy)
    for (size_t b = 0; b < nnb; ++b) {
        int32_t const & nb = nbs[b];
        int64_t nbt = nbtarget[b];
        if (nbt < 0) {
            // This neighbor is disabled.  Check for collisions with the neighbor
            // in its fixed theta / phi location.
            collide = hw->collide_xy_thetaphi(
                loc,
                tpos,
                nb,
                hw->loc_theta_pos.at(nb),
                hw->loc_phi_pos.at(nb)
            );
        } else {
            // Neighbor is working, check for collisions with the neighbor in
            // its currently assigned position.
            auto npos = target_xy.at(nbt);
            collide = hw->collide_xy(loc, tpos, nb, npos);
        }
        // Remove these lines if switching back to threading.
        if (collide) {
            if (extra_log) {
                logmsg.str("");
                logmsg << "ok_to_assign: tile " << tile << ", loc "
                    << loc << ", target " << target
                    << " would collide with target " << nbt;
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
            return false;
        }
    }

    // Restore these lines if switching back to threading.
    // if (collide) {
    //     return false;
    // }

    // Would this assignment hit a GFA or petal edge?

    collide = hw->collide_xy_edges(loc, tpos);
    if (collide) {
        if (extra_log) {
            logmsg.str("");
            logmsg << "ok_to_assign: tile " << tile << ", loc "
                << loc << ", target " << target
                << " would collide with GFA or Petal Boundary ";
            logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
        }
        return false;
    }

    // All good!
    return true;
}


void fba::Assignment::assign_tileloc(fba::Hardware const * hw,
    fba::Targets * tgs, int32_t tile, int32_t loc, int64_t target,
    uint8_t type) {

    fba::Logger & logger = fba::Logger::get();
    bool extra_log = logger.extra_debug();
    std::ostringstream logmsg;

    if (target < 0) {
        logmsg.str("");
        logmsg << "cannot assign negative target ID to tile "
            << tile << ", loc " << loc << ".  Did you mean to unassign?";
        logger.warning(logmsg.str().c_str());
        return;
    }

    auto & ftarg = loc_target[tile];

    if (ftarg.count(loc) > 0) {
        int64_t cur = ftarg.at(loc);
        if (cur >= 0) {
            logmsg.str("");
            logmsg << "tile " << tile << ", loc " << loc
                << " already assigned to target " << cur
                << " cannot assign " << target;
            logger.warning(logmsg.str().c_str());
            return;
        }
    }

    auto & tfiber = target_loc[target];

    if (tfiber.count(tile) > 0) {
        logmsg.str("");
        logmsg << "target " << target << " already assigned on tile " << tile;
        logger.warning(logmsg.str().c_str());
        return;
    }

    auto & tgobj = tgs->data.at(target);

    if ( ! tgobj.is_type(type)) {
        logmsg.str("");
        logmsg << "target " << target << " not of type "
            << (int)type;
        logger.error(logmsg.str().c_str());
        throw std::runtime_error(logmsg.str().c_str());
    }

    ftarg[loc] = target;
    target_loc[target][tile] = loc;

    int32_t petal = hw->loc_petal.at(loc);
    int32_t slitblock = hw->loc_slitblock.at(loc);

    // Objects can be more than one type (e.g. standards and science).  When
    // incrementing the counts of object types per tile and petal, we want
    // to update the counts for valid types of this object.
    static const std::vector <uint8_t> target_types = {
        TARGET_TYPE_SCIENCE, TARGET_TYPE_STANDARD,
        TARGET_TYPE_SKY, TARGET_TYPE_SUPPSKY, TARGET_TYPE_SAFE};
    for (auto const & tt : target_types) {
        if (tgobj.is_type(tt)) {
            nassign_tile.at(tt).at(tile)++;
            nassign_petal.at(tt).at(tile).at(petal)++;
            nassign_slitblock.at(tt).at(tile).at(petal).at(slitblock)++;
            if (extra_log) {
                logmsg.str("");
                logmsg << "assign_tileloc: tile " << tile << ", loc "
                    << loc << ", target " << target << ", type "
                    << (int)tt << " N_tile now = "
                    << nassign_tile.at(tt).at(tile)
                    << " N_petal now = "
                    << nassign_petal.at(tt).at(tile).at(petal)
                    << " N_slitblock now = "
                    << nassign_slitblock.at(tt).at(tile).at(petal).at(slitblock);
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
        }
    }
    tgobj.obsremain--;

    return;
}


void fba::Assignment::unassign_tileloc(fba::Hardware const * hw,
    fba::Targets * tgs, int32_t tile, int32_t loc, uint8_t type) {

    fba::Logger & logger = fba::Logger::get();
    bool extra_log = logger.extra_debug();
    std::ostringstream logmsg;

    if (loc_target.count(tile) == 0) {
        logmsg.str("");
        logmsg << "tile " << tile
            << " has no locations assigned.  Ignoring unassign";
        logger.warning(logmsg.str().c_str());
        return;
    }
    auto & ftarg = loc_target.at(tile);

    if (ftarg.count(loc) == 0) {
        logmsg.str("");
        logmsg << "tile " << tile << ", loc " << loc
            << " already unassigned";
        logger.warning(logmsg.str().c_str());
        return;
    }

    int64_t target = ftarg.at(loc);
    if (target < 0) {
        logmsg.str("");
        logmsg << "tile " << tile << ", loc " << loc
            << " already unassigned";
        logger.warning(logmsg.str().c_str());
        return;
    }

    auto & tgobj = tgs->data.at(target);

    if ( ! tgobj.is_type(type)) {
        logmsg.str("");
        logmsg << "current target " << target << " not of type "
            << (int)type << " requested in unassign of tile " << tile
            << ", loc " << loc;
        logger.error(logmsg.str().c_str());
        throw std::runtime_error(logmsg.str().c_str());
    }

    int32_t petal = hw->loc_petal.at(loc);
    int32_t slitblock = hw->loc_slitblock.at(loc);

    // Objects can be more than one type (e.g. standards and science).  When
    // incrementing the counts of object types per tile and petal, we want
    // to update the counts for valid types of this object.
    static const std::vector <uint8_t> target_types = {
        TARGET_TYPE_SCIENCE, TARGET_TYPE_STANDARD,
        TARGET_TYPE_SKY, TARGET_TYPE_SUPPSKY, TARGET_TYPE_SAFE};
    for (auto const & tt : target_types) {
        if (tgobj.is_type(tt)) {
            nassign_tile.at(tt).at(tile)--;
            nassign_petal.at(tt).at(tile).at(petal)--;
            nassign_slitblock.at(tt).at(tile).at(petal).at(slitblock)--;
            if (extra_log) {
                logmsg.str("");
                logmsg << "unassign_tileloc: tile " << tile << ", loc "
                    << loc << ", target " << target << ", type "
                    << (int)tt << " N_tile now = "
                    << nassign_tile.at(tt).at(tile)
                    << " N_petal now = "
                    << nassign_petal.at(tt).at(tile).at(petal)
                    << " N_slitblock now = "
                    << nassign_slitblock.at(tt).at(tile).at(petal).at(slitblock);
                logger.debug_tfg(tile, loc, target, logmsg.str().c_str());
            }
        }
    }
    tgobj.obsremain++;

    target_loc[target].erase(tile);
    ftarg.erase(loc);

    return;
}


void fba::Assignment::targets_to_project(
    fba::Targets const * tgs,
    std::map <int32_t, std::vector <int64_t> > const & tgsavail,
    std::vector <int32_t> const & locs,
    std::vector <int64_t> & tgids,
    std::vector <double> & tgra,
    std::vector <double> & tgdec) const {
    // This function computes the target IDs that need to be projected
    // for a given set of locations on a tile.

    std::set <int64_t> seen;
    tgids.clear();
    tgra.clear();
    tgdec.clear();

    for (auto const & lid : locs) {
        // Project this location's targets
        if (tgsavail.count(lid) > 0) {
            // The available targets for this location.
            auto const & avail = tgsavail.at(lid);
            for (auto const & id : avail) {
                if (seen.count(id) == 0) {
                    // This target has not yet been processed.
                    auto const & tg = tgs->data.at(id);
                    tgids.push_back(id);
                    tgra.push_back(tg.ra);
                    tgdec.push_back(tg.dec);
                    seen.insert(id);
                }
            }
        }
    }

    return;
}


void fba::Assignment::project_targets(fba::Hardware const * hw,
        fba::Targets const * tgs, fba::TargetsAvailable const * tgsavail,
        int32_t tile_id, double tile_ra, double tile_dec, double tile_theta,
        std::map <int64_t, std::pair <double, double> > & target_xy) const {
    // This function computes the projection of all available targets
    // into focalplane coordinates for one tile.  We do this once and then
    // use these positions multiple times.

    // Clear the output
    target_xy.clear();

    if (tgsavail->data.count(tile_id) == 0) {
        // This tile has no locations with available targets.
        return;
    }

    // Reference to the available targets for this tile.
    auto const & tfavail = tgsavail->data.at(tile_id);

    if (tfavail.size() == 0) {
        // There are no locations on this tile with available targets.
        return;
    }

    // List of locations we are projecting- anything with targets available.
    std::vector <int32_t> lids;
    for (auto const & f : hw->locations) {
        if ((tfavail.count(f) != 0) && (tfavail.at(f).size() > 0)) {
            lids.push_back(f);
        }
    }

    // Vectors of all targets to compute
    std::vector <int64_t> tgids;
    std::vector <double> tgra;
    std::vector <double> tgdec;

    targets_to_project(tgs, tfavail, lids, tgids, tgra, tgdec);

    // Now thread over the targets to compute.  Note we are always working in the
    // curved focal surface for assignment.
    bool use_CS5 = false;

    std::vector <std::pair <double, double> > xy;
    hw->radec2xy_multi(tile_ra, tile_dec, tile_theta, tgra, tgdec, xy, use_CS5);

    for (size_t t = 0; t < tgids.size(); ++t) {
        target_xy[tgids[t]] = xy[t];
    }

    return;
}



fba::Hardware::pshr fba::Assignment::hardware() const {
    return hw_;
}


fba::Targets::pshr fba::Assignment::targets() const {
    return tgs_;
}


fba::Tiles::pshr fba::Assignment::tiles() const {
    return tiles_;
}


fba::TargetsAvailable::pshr fba::Assignment::targets_avail() const {
    return tgsavail_;
}


fba::LocationsAvailable::pshr fba::Assignment::locations_avail() const {
    return locavail_;
}


// void fba::Assignment::dump_fits(std::string const & prefix) const {
//     fba::Timer tm;
//     tm.start();
//
//     fba::GlobalTimers & gtm = fba::GlobalTimers::get();
//     std::ostringstream gtmname;
//
//     fba::Logger & logger = fba::Logger::get();
//     std::ostringstream logmsg;
//
//     gtmname.str("");
//     gtmname << "Assignment dump_fits";
//     gtm.start(gtmname.str());
//
//     fba::Tiles::pshr tiles = tgsavail_->tiles();
//     size_t ntile = tiles->id.size();
//
//     auto const * ptiles = tiles.get();
//     auto const * phw = hw_.get();
//     auto const * ptgs = tgs_.get();
//     auto const * ptgsavail = tgsavail_.get();
//
//     #pragma omp parallel for schedule(dynamic) default(none) shared(ntile, ptiles, phw, ptgs, ptgsavail, tile_target_xy, logmsg, logger)
//     for (size_t t = 0; t < ntile; ++t) {
//         int32_t tile_id = ptiles->id[t];
//         double tile_ra = ptiles->ra[t];
//         double tile_dec = ptiles->dec[t];
//         char tile_file[1024];
//         int ret = snprintf(tile_file, 1024, "%s_%06d.fits", prefix, tile_id);
//
//
//
//     }
//
//     gtm.stop(gtmname.str());
//
//     gtmname.str("");
//     gtmname << "Assignment dump_fits";
//     gtm.stop(gtmname.str());
//
//     logmsg.str("");
//     logmsg << "Assignment dump_fits";
//     tm.stop();
//     tm.report(logmsg.str().c_str());
//
//     return;
// }
