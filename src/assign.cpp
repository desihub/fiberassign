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
                            fba::FibersAvailable::pshr favail) {
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
    favail_ = favail;

    // Get the hardware and tile configuration
    tiles_ = tgsavail_->tiles();
    hw_ = tgsavail_->hardware();

    // Initialize assignment counts

    std::vector <uint8_t> tgtypes;
    tgtypes.push_back(TARGET_TYPE_SCIENCE);
    tgtypes.push_back(TARGET_TYPE_STANDARD);
    tgtypes.push_back(TARGET_TYPE_SKY);
    tgtypes.push_back(TARGET_TYPE_SAFE);

    tile_target_xy.clear();

    size_t ntile = tiles_->id.size();
    for (auto const & tp : tgtypes) {
        nassign_tile[tp].clear();
        nassign_petal[tp].clear();
        for (size_t t = 0; t < ntile; ++t) {
            int32_t tile_id = tiles_->id[t];
            nassign_tile[tp][tile_id] = 0;
            nassign_petal[tp][tile_id].clear();
            for (int32_t p = 0; p < hw_->npetal; ++p) {
                nassign_petal[tp][tile_id][p] = 0;
            }
            tile_target_xy[tile_id].clear();
        }
    }

    fiber_target.clear();
    target_fiber.clear();

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
        std::map <int64_t, std::pair <double, double> > local_xy;
        project_targets(phw, ptgs, ptgsavail, tile_id, tile_ra, tile_dec,
                        local_xy);
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
    for (auto const & it : fiber_target) {
        ret.push_back(it.first);
    }
    return ret;
}


std::map <int32_t, int64_t> const & fba::Assignment::tile_fiber_target(int32_t tile) const {
    return fiber_target.at(tile);
}


void fba::Assignment::parse_tile_range(int32_t start_tile, int32_t stop_tile,
                                       int32_t & tstart, int32_t & tstop) {
    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    int32_t ntile = tiles_->id.size();

    tstart = -1;
    tstop = -1;
    if (start_tile >= 0) {
        for (int32_t t = 0; t < ntile; ++t) {
            if (tiles_->id[t] == start_tile) {
                tstart = t;
            }
        }
    }
    if (tstart < 0) {
        // was not found
        tstart = 0;
        if (start_tile >= 0) {
            logmsg.str("");
            logmsg << "requested start tile " << start_tile
                << " was not found in tile list.  Starting at tile "
                << tiles_->id[tstart] << " instead";
            logger.warning(logmsg.str().c_str());
        }
    }
    if (stop_tile >= 0) {
        for (int32_t t = 0; t < ntile; ++t) {
            if (tiles_->id[t] == stop_tile) {
                tstop = t;
            }
        }
    }
    if (tstop < 0) {
        // was not found
        tstop = ntile - 1;
        if (stop_tile >= 0) {
            logmsg.str("");
            logmsg << "requested stop tile " << stop_tile
                << " was not found in tile list.  Stopping at tile "
                << tiles_->id[tstop] << " instead";
            logger.warning(logmsg.str().c_str());
        }
    }

    return;
}


void fba::Assignment::assign_unused(uint8_t tgtype, int32_t max_per_petal,
                                    std::string const & pos_type,
                                    int32_t start_tile, int32_t stop_tile) {
    fba::Timer tm;
    tm.start();

    fba::GlobalTimers & gtm = fba::GlobalTimers::get();
    std::ostringstream gtmname;

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    std::string tgstr = fba::target_string(tgtype);

    // Select fiber IDs based on positioner type
    auto fiber_id = hw_->device_fibers(pos_type);
    logmsg.str("");
    logmsg << "assign unused " << tgstr << ":  considering "
        << fiber_id.size()
        << " fibers of positioner type \"" << pos_type << "\"";
    logger.info(logmsg.str().c_str());

    gtmname.str("");
    gtmname << "unused " << tgstr << ": total";
    gtm.start(gtmname.str());

    if (max_per_petal < 0) {
        // A negative value indicates that there is no limit.
        max_per_petal = 2147483647;
    }

    // Determine our range of tiles
    int32_t tstart;
    int32_t tstop;
    parse_tile_range(start_tile, stop_tile, tstart, tstop);

    logmsg.str("");
    logmsg << "assign unused " << tgstr << ":  working on tiles "
        << start_tile << " (index " << tstart << ") to "
        << stop_tile << " (index " << tstop << ")";
    logger.info(logmsg.str().c_str());

    auto const & tavail = tgsavail_->data;
    auto const & tgsdata = tgs_->data;

    fba::weight_compare fiber_comp;

    // This container stores the position in focalplane
    // coordinates of the available targets on a tile.
    //std::map <int64_t, std::pair <double, double> > target_xy;

    // This is the specific list of available targets for a single
    // tile / fiber, after selecting by target type, etc.
    std::vector <int64_t> targets_avail;

    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];
        double tile_ra = tiles_->ra[t];
        double tile_dec = tiles_->dec[t];

        logmsg.str("");
        logmsg << "assign unused " << tgstr << ": working on tile " << tile_id
            << " at RA/DEC = " << tile_ra << " / " << tile_dec;
        logger.debug(logmsg.str().c_str());

        if ((tavail.count(tile_id) == 0)
            || (tavail.at(tile_id).size() == 0)) {
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

        // Available targets for this tile.
        auto const & tfavail = tavail.at(tile_id);
        auto const & target_xy = tile_target_xy.at(tile_id);

        // The order in which the fibers should be assigned, based
        // on the maximum priority of each fiber's available targets of the
        // specified type.
        std::vector <weight_index> fiber_priority;

        gtmname.str("");
        gtmname << "unused " << tgstr << ": fiber order";
        gtm.start(gtmname.str());

        for (auto const & fid : fiber_id) {
            if (tfavail.count(fid) == 0) {
                // No targets available for this fiber
                continue;
            }
            // Targets available for this fiber.  These are already
            // sorted by priority + subpriority.
            auto const & avail = tfavail.at(fid);

            if (avail.size() == 0) {
                // No targets available for this fiber.
                continue;
            }

            // The list of available targets is already sorted by
            // total priority.  Just grab the first target that is the
            // correct type and has remaining observations.
            for (auto const & tgid : avail) {
                // Reference to the Target object with this ID
                auto const & tg = tgsdata.at(tgid);

                if ( ! tg.is_type(tgtype)) {
                    // This is not the correct target type.
                    continue;
                }

                if ((tgtype == TARGET_TYPE_SCIENCE) && (tg.obs_remain <= 0)) {
                    // Done observing science observations for this target
                    continue;
                }

                double totpriority = static_cast <double> (tg.priority)
                    + tg.subpriority;

                fiber_priority.push_back(std::make_pair(totpriority, fid));
                break;
            }
        }

        std::stable_sort(fiber_priority.begin(), fiber_priority.end(),
                         fiber_comp);
        gtm.stop(gtmname.str());

        gtmname.str("");
        gtmname << "unused " << tgstr << ": assign fibers";
        gtm.start(gtmname.str());

        for (auto const & fpr : fiber_priority) {
            // The fiber ID.
            int32_t fid = fpr.second;
            // The petal
            int32_t p = hw_->fiber_petal[fid];

            if (nassign_petal[tgtype][tile_id][p] >= max_per_petal) {
                // Already have enough objects on this petal
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "assign unused " << tgstr << ": tile " << tile_id
                        << ", petal " << p << " fiber " << fid << " has "
                        << nassign_petal[tgtype][tile_id][p]
                        << " (>= " << max_per_petal << ")";
                    logger.debug_tfg(tile_id, fid, -1, logmsg.str().c_str());
                }
                continue;
            }

            if (fiber_target[tile_id].count(fid) > 0) {
                // Fiber is currently assigned, skip it.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "assign unused " << tgstr << ": tile " << tile_id
                        << ", petal " << p << " fiber " << fid
                        << " is already assigned";
                    logger.debug_tfg(tile_id, fid, -1, logmsg.str().c_str());
                }
                continue;
            }

            // Examine available targets for this fiber.  These are already
            // sorted in priority / subpriority order.  Attempt to assign
            // any standards we can.

            if (tfavail.count(fid) == 0) {
                // No targets available for this fiber
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "assign unused " << tgstr << ": tile " << tile_id
                        << ", petal " << p << " fiber " << fid
                        << " has no available targets";
                    logger.debug_tfg(tile_id, fid, -1, logmsg.str().c_str());
                }
                continue;
            }

            // All targets available for this fiber.
            auto const & avail = tfavail.at(fid);

            if (avail.size() == 0) {
                // No targets available for this fiber.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "assign unused " << tgstr << ": tile " << tile_id
                        << ", petal " << p << " fiber " << fid
                        << " has no available targets";
                    logger.debug_tfg(tile_id, fid, -1, logmsg.str().c_str());
                }
                continue;
            }

            // Find available targets of the desired type

            targets_avail.clear();
            for (auto const & tgid : avail) {
                if (tgsdata.count(tgid) == 0) {
                    std::cerr << "available target " << tgid << " does not exist in target list!" << std::endl;
                    throw std::runtime_error("missing target");
                }
                // Reference to the Target object with this ID
                auto const & tg = tgsdata.at(tgid);
                if ( ! tg.is_type(tgtype)) {
                    // This is not the desired type
                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign unused " << tgstr << ": tile " << tile_id
                            << ", petal " << p << " fiber " << fid
                            << " available target " << tgid
                            << " is wrong type (" << (int)tg.type << ")";
                        logger.debug_tfg(tile_id, fid, tgid,
                                         logmsg.str().c_str());
                    }
                    continue;
                }
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "assign unused " << tgstr << ": tile " << tile_id
                        << ", petal " << p << " fiber " << fid
                        << " available target " << tgid
                        << ", priority " << tg.priority << ", subpriority "
                        << tg.subpriority;
                    logger.debug_tfg(tile_id, fid, tgid,
                                     logmsg.str().c_str());
                }
                targets_avail.push_back(tgid);
            }

            // Assign the best object that is not already assigned.

            if (targets_avail.size() > 0) {
                int64_t target = find_best(hw_.get(), tgs_.get(), tile_id,
                                           fid, tgtype, target_xy, targets_avail);
                if (target >= 0) {
                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign unused " << tgstr << ": tile " << tile_id
                            << ", petal " << p << " fiber " << fid
                            << " found best object " << target;
                        logger.debug_tfg(tile_id, fid, target,
                                         logmsg.str().c_str());
                    }
                    assign_tilefiber(hw_.get(), tgs_.get(), tile_id,
                                     fid, target, tgtype);
                }
            } else {
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "assign unused " << tgstr << ": tile " << tile_id
                        << ", petal " << p << " fiber " << fid
                        << " no available targets of correct type";
                    logger.debug_tfg(tile_id, fid, -1,
                                     logmsg.str().c_str());
                }
            }
        }

        gtm.stop(gtmname.str());
    }

    gtmname.str("");
    gtmname << "unused " << tgstr << ": total";
    gtm.stop(gtmname.str());

    logmsg.str("");
    if (max_per_petal > 1000000) {
        logmsg << "Assign " << tgstr << " targets to unused fibers";
    } else {
        logmsg << "Assign up to " << max_per_petal << " " << tgstr
            << " targets to unused fibers";
    }
    tm.stop();
    tm.report(logmsg.str().c_str());

    return;
}


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

    // Select fiber IDs based on positioner type
    auto fiber_id = hw_->device_fibers("POS");

    // Determine our range of tiles
    int32_t tstart;
    int32_t tstop;
    parse_tile_range(start_tile, stop_tile, tstart, tstop);

    logmsg.str("");
    logmsg << "redist:  working on tiles "
        << start_tile << " (index " << tstart << ") to "
        << stop_tile << " (index " << tstop << ")";
    logger.info(logmsg.str().c_str());

    // This is going to be used to track whether a tile / fiber combination has
    // already been considered.
    std::map <int32_t, std::map <int32_t, bool> > done;
    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];
        done[tile_id].clear();
    }

    std::vector <int64_t> targets_avail;

    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];

        if (nassign_tile[TARGET_TYPE_SCIENCE][tile_id] == 0) {
            // Skip tiles that are fully unassigned, since this implies that the
            // first pass of the assignment had no available targets to use.
            continue;
        }

        logmsg.str("");
        logmsg << "redist: working on tile " << tile_id
            << " with " << nassign_tile[TARGET_TYPE_SCIENCE][tile_id]
            << " science targets currently assigned";
        logger.debug(logmsg.str().c_str());

        auto const & tfavail = tgsavail_->data.at(tile_id);
        auto const & target_xy = tile_target_xy.at(tile_id);

        for (auto const & fid : fiber_id) {
            if (fiber_target[tile_id].count(fid) == 0) {
                // This tile / fiber combination is unassigned.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "redist: tile " << tile_id << ", fiber "
                        << fid << " unassigned- skipping";
                    logger.debug_tfg(tile_id, fid, -1, logmsg.str().c_str());
                }
                continue;
            }

            if ((done[tile_id].count(fid) > 0) && done[tile_id].at(fid)) {
                // Already considered or swapped this fiber.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "redist: tile " << tile_id << ", fiber "
                        << fid << " already considered- skipping";
                    logger.debug_tfg(tile_id, fid, -1, logmsg.str().c_str());
                }
                continue;
            }

            // Get the current target.
            int64_t tgid = fiber_target[tile_id].at(fid);
            auto const & tg = tgs_->data.at(tgid);
            if ( ! tg.is_science()) {
                // Only consider science targets.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "redist: tile " << tile_id << ", fiber "
                        << fid << ", target " << tgid
                        << " not a science target- skipping";
                    logger.debug_tfg(tile_id, fid, tgid, logmsg.str().c_str());
                }
                continue;
            }

            // Find any better assignment.
            int32_t best_tile;
            int32_t best_fiber;

            reassign_science_target(tstart, tstop, tile_id, fid, tgid, true,
                                    done, best_tile, best_fiber);

            // Mark this current tile / fiber as done
            done[tile_id][fid] = true;

            if ((best_tile != tile_id) || (best_fiber != fid)) {
                // We have a better possible assignment- change it.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "redist: tile " << tile_id << ", fiber "
                        << fid << ", target " << tgid
                        << " swapping to T/F " << best_tile << ","
                        << best_fiber;
                    logger.debug_tfg(tile_id, fid, tgid, logmsg.str().c_str());
                }
                unassign_tilefiber(hw_.get(), tgs_.get(), tile_id, fid,
                                   TARGET_TYPE_SCIENCE);
                assign_tilefiber(hw_.get(), tgs_.get(), best_tile,
                                 best_fiber, tgid, TARGET_TYPE_SCIENCE);
                // Mark the new assignment as done.
                done[best_tile][best_fiber] = true;

                // All targets available for this fiber.
                auto const & avail = tfavail.at(fid);

                if (avail.size() > 0) {
                    // Find available targets of the desired type
                    targets_avail.clear();
                    for (auto const & tgid : avail) {
                        // Reference to the Target object with this ID
                        auto const & tg = tgs_->data.at(tgid);
                        if (tg.is_science()) {
                            targets_avail.push_back(tgid);
                        }
                    }
                    // Assign the best object that is not already assigned.
                    if (targets_avail.size() > 0) {
                        int64_t newtarget = find_best(hw_.get(), tgs_.get(),
                                                      tile_id, fid,
                                                      TARGET_TYPE_SCIENCE, target_xy,
                                                      targets_avail);
                        if (newtarget >= 0) {
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "redist: tile " << tile_id
                                    << " fiber " << fid
                                    << " reassigned to " << newtarget;
                                logger.debug_tfg(tile_id, fid, newtarget,
                                                 logmsg.str().c_str());
                            }
                            assign_tilefiber(hw_.get(), tgs_.get(), tile_id,
                                             fid, newtarget,
                                             TARGET_TYPE_SCIENCE);
                        } else {
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "redist: tile " << tile_id
                                    << " fiber " << fid
                                    << " no science targets available ";
                                logger.debug_tfg(tile_id, fid, newtarget,
                                                 logmsg.str().c_str());
                            }
                        }
                    }
                }
            } else {
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "redist: tile " << tile_id << ", fiber "
                        << fid << ", target " << tgid
                        << " keeping assignment";
                    logger.debug_tfg(tile_id, fid, tgid, logmsg.str().c_str());
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


// Helper types for sorting targets by subpriority

typedef std::pair <double, int64_t> weight_target;

struct target_compare {
    bool operator() (weight_target const & lhs,
                     weight_target const & rhs) const {
        return lhs.first > rhs.first;
    }
};

typedef std::pair <double, int32_t> weight_fiber;

struct fiber_compare {
    bool operator() (weight_fiber const & lhs,
                     weight_fiber const & rhs) const {
        return lhs.first < rhs.first;
    }
};


void fba::Assignment::assign_force(uint8_t tgtype, int32_t required_per_petal,
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

    int32_t npetal = hw_->npetal;
    auto fiber_id = hw_->device_fibers("POS");

    // Determine our range of tiles
    int32_t tstart;
    int32_t tstop;
    parse_tile_range(start_tile, stop_tile, tstart, tstop);

    logmsg.str("");
    logmsg << "assign force " << tgstr << ":  working on tiles "
        << start_tile << " (index " << tstart << ") to "
        << stop_tile << " (index " << tstop << ")";
    logger.info(logmsg.str().c_str());

    auto const & tgsdata = tgs_->data;

    // This is going to be used to track whether a tile / fiber combination has
    // already been considered.
    std::map <int32_t, std::map <int32_t, bool> > done;
    for (int32_t t = tstart; t <= tstop; ++t) {
        int32_t tile_id = tiles_->id[t];
        done[tile_id].clear();
    }

    // The target classes
    auto const & tgclasses = tgs_->science_classes;

    // This is the specific list of available targets for a single
    // tile / fiber, after selecting by target type, etc.
    std::vector <int64_t> targets_avail;

    // vector of priority-weighted objects on one petal
    std::vector <weight_target> petal_obj;

    // The fibers on one petal, of each priority class, that can reach
    // each object.
    std::vector <weight_fiber> can_replace;

    target_compare tcompare;
    fiber_compare fcompare;

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

        // Available targets for this tile.
        auto const & tfavail = tgsavail_->data.at(tile_id);
        auto const & target_xy = tile_target_xy.at(tile_id);

        for (int32_t p = 0; p < npetal; ++p) {
            if (extra_log) {
                logmsg.str("");
                logmsg << "assign force " << tgstr << ": tile " << tile_id
                    << ", petal " << p << " has "
                    << nassign_petal[tgtype][tile_id][p]
                    << " (require " << required_per_petal << ")";
                logger.debug_tfg(tile_id, -1, -1, logmsg.str().c_str());
            }
            if (nassign_petal[tgtype][tile_id][p] >= required_per_petal) {
                // Already have enough objects on this petal
                continue;
            }

            // Targets of desired type on this petal.

            gtmname.str("");
            gtmname << "force " << tgstr << ": per petal possible objects";
            gtm.start(gtmname.str());

            // This object tracks the available standards/sky on this petal
            // which can replace targets of the current target class.  These
            // are sorted by total priority of the standard/sky.
            petal_obj.clear();

            // Get the fibers for this petal which are science positioners.
            std::vector <int32_t> petal_fibers;
            for (auto const & fid : fiber_id) {
                if (hw_->fiber_petal.at(fid) == p) {
                    petal_fibers.push_back(fid);
                }
            }

            for (auto const & fid : petal_fibers) {
                if (tfavail.count(fid) == 0) {
                    // No targets available for this fiber
                    continue;
                }

                // Targets available for this fiber.
                auto const & avail = tfavail.at(fid);

                if (avail.size() == 0) {
                    continue;
                }

                for (auto const & tgid : avail) {
                    // Reference to the Target object with this ID
                    auto & tg = tgsdata.at(tgid);
                    if ((tg.type & tgtype) == 0) {
                        // This is an object of the wrong type.
                        continue;
                    }

                    double av_weight = static_cast <double> (tg.priority)
                        + tg.subpriority;
                    petal_obj.push_back(std::make_pair(av_weight, tgid));

                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign force " << tgstr << ": tile "
                            << tile_id << ", petal " << p << ", fiber "
                            << fid << ", found object " << tgid
                            << " with weight " << av_weight;
                        logger.debug_tfg(tile_id, fid, tgid,
                                         logmsg.str().c_str());
                    }
                }
            }

            if (extra_log) {
                logmsg.str("");
                logmsg << "assign force " << tgstr << ": tile "
                    << tile_id << ", petal " << p
                    << " has " << petal_obj.size() << " available objects";
                logger.debug_tfg(tile_id, -1, -1, logmsg.str().c_str());
            }
            std::stable_sort(petal_obj.begin(), petal_obj.end(), tcompare);

            gtm.stop(gtmname.str());

            // Go through all target priority classes in reverse order.

            for (auto const & tc : tgclasses) {

                // Go through new objects in priority order and try to bump
                // science targets of the given class.
                for (auto const & ps : petal_obj) {
                    // target ID for this object
                    int64_t id = ps.second;

                    // Find all fibers on this petal which can reach
                    // this target and are occupied by a science target of
                    // the correct target class.  Sort these by increasing
                    // subpriority.
                    can_replace.clear();
                    for (auto const & tf : favail_->data[id]) {
                        int32_t fid = tf.second;
                        if (tf.first != tile_id) {
                            // Not this tile
                            continue;
                        }
                        if (hw_->fiber_petal[fid] != p) {
                            // Not this petal
                            continue;
                        }
                        // Is this fiber currently assigned to a target of the
                        // correct target class?
                        if (fiber_target[tile_id].count(fid) == 0) {
                            // Fiber not assigned
                            continue;
                        }
                        auto & cur = tgsdata.at(fiber_target[tile_id].at(fid));
                        if (cur.priority != tc) {
                            // Wrong target class
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "assign force " << tgstr << ": tile "
                                    << tile_id << ", petal " << p
                                    << ", class " << tc << ", object " << id
                                    << ", subpriority " << ps.first
                                    << ", available fiber " << fid
                                    << " at target " << cur.id
                                    << " is wrong class (" << cur.priority
                                    << ")";
                                logger.debug_tfg(tile_id, fid, id,
                                                 logmsg.str().c_str());
                            }
                            continue;
                        }

                        // Special case:  if the target is both a science
                        // target AND a standard, then we never bump it.  The
                        // reason is the following:  if we are bumping science
                        // for standards, then this is self defeating.  If we
                        // are bumping science targets for sky, then this
                        // would change previous per-petal standards counts.
                        // So we just never do it.

                        if (cur.is_science() && cur.is_standard()) {
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "assign force " << tgstr << ": tile "
                                    << tile_id << ", petal " << p
                                    << ", class " << tc << ", object " << id
                                    << ", subpriority " << ps.first
                                    << ", available fiber " << fid
                                    << " at science target " << cur.id
                                    << " is also a standard- skipping";
                                logger.debug_tfg(tile_id, fid, id,
                                                 logmsg.str().c_str());
                            }
                            continue;
                        }
                        can_replace.push_back(
                            std::make_pair(cur.subpriority, fid));
                    }
                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign force " << tgstr << ": tile "
                            << tile_id << ", petal " << p
                            << ", class " << tc << ", object " << id
                            << ", subpriority " << ps.first
                            << " has " << can_replace.size()
                            << " potential targets to bump";
                        logger.debug_tfg(tile_id, -1, id,
                                         logmsg.str().c_str());
                    }

                    std::stable_sort(can_replace.begin(), can_replace.end(),
                                     fcompare);

                    for (auto const & av : can_replace) {
                        int32_t fid = av.second;
                        if ((done[tile_id].count(fid) > 0)
                            && done[tile_id].at(fid)) {
                            // Already considered this tile/fiber
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "assign force " << tgstr << ": tile "
                                    << tile_id << ", petal " << p
                                    << ", class " << tc << ", object " << id
                                    << ", subpriority " << ps.first
                                    << ", available fiber " << fid
                                    << " already bumped";
                                logger.debug_tfg(tile_id, fid, id,
                                                 logmsg.str().c_str());
                            }
                            continue;
                        }

                        int64_t curtg = fiber_target[tile_id].at(fid);

                        if (ok_to_assign(hw_.get(), tile_id, fid,
                                         id, target_xy)) {
                            // We can swap this fiber.
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "assign force " << tgstr << ": tile "
                                    << tile_id << ", petal " << p
                                    << ", class " << tc << ", object " << id
                                    << ", subpriority " << ps.first
                                    << ", available fiber " << fid
                                    << " bumping science target " << curtg;
                                // Log for target doing the bumping
                                logger.debug_tfg(tile_id, fid, id,
                                                 logmsg.str().c_str());
                                // Also log for target getting bumped
                                logger.debug_tfg(tile_id, fid, curtg,
                                                 logmsg.str().c_str());
                            }
                            // Attempt to re-assign this science target to
                            // a later tile / fiber.
                            int32_t best_tile;
                            int32_t best_fiber;
                            reassign_science_target(tstart, tstop, tile_id,
                                                    fid, curtg, false,
                                                    done, best_tile,
                                                    best_fiber);

                            // Assign object.
                            unassign_tilefiber(hw_.get(), tgs_.get(),
                                               tile_id, fid,
                                               TARGET_TYPE_SCIENCE);
                            assign_tilefiber(hw_.get(), tgs_.get(),
                                             tile_id, fid, id, tgtype);

                            // Mark current tile / fiber as done
                            done[tile_id][fid] = true;

                            // Move science target to new location if possible
                            if ((best_tile != tile_id)
                                && (best_fiber != fid)) {
                                assign_tilefiber(hw_.get(), tgs_.get(),
                                                 best_tile, best_fiber,
                                                 curtg, TARGET_TYPE_SCIENCE);
                                // Mark this new location as done
                                done[best_tile][best_fiber] = true;
                                if (extra_log) {
                                    logmsg.str("");
                                    logmsg << "assign force " << tgstr
                                        << ": tile " << tile_id
                                        << ", petal " << p
                                        << ", class " << tc << ", object " << id
                                        << ", subpriority " << ps.first
                                        << ", available fiber " << fid
                                        << " bumped target " << curtg
                                        << " reassigned to " << best_tile
                                        << "," << best_fiber;
                                    // Log for target doing the bumping
                                    logger.debug_tfg(tile_id, fid, id,
                                                     logmsg.str().c_str());
                                    // Also log for target getting bumped
                                    logger.debug_tfg(tile_id, fid, curtg,
                                                     logmsg.str().c_str());
                                }
                            }
                            break;
                        } else {
                            if (extra_log) {
                                logmsg.str("");
                                logmsg << "assign force " << tgstr
                                    << ": tile " << tile_id
                                    << ", petal " << p
                                    << ", class " << tc << ", object " << id
                                    << ", subpriority " << ps.first
                                    << ", available fiber " << fid
                                    << " not ok to assign";
                                logger.debug_tfg(tile_id, fid, id,
                                                 logmsg.str().c_str());
                            }
                        }
                    }

                    if (extra_log) {
                        logmsg.str("");
                        logmsg << "assign force " << tgstr
                            << ": tile " << tile_id
                            << ", petal " << p << " now has "
                            << nassign_petal[tgtype][tile_id][p]
                            << " (require " << required_per_petal << ")";
                        logger.debug_tfg(tile_id, -1, -1, logmsg.str().c_str());
                    }
                    if (nassign_petal[tgtype][tile_id][p]
                        >= required_per_petal) {
                        // Have enough objects on this petal
                        break;
                    }
                }
                if (nassign_petal[tgtype][tile_id][p]
                    >= required_per_petal) {
                    // Have enough objects on this petal
                    break;
                }
            }

            if (nassign_petal[tgtype][tile_id][p]
                < required_per_petal) {
                logmsg.str("");
                logmsg << "assign force " << tgstr << ": tile " << tile_id
                    << ", petal " << p << " could only assign "
                    << nassign_petal[tgtype][tile_id][p]
                    << " (require " << required_per_petal
                    << ").  Insufficient number of objects or too many collisions";
                logger.warning(logmsg.str().c_str());
            }
        }
    }

    gtmname.str("");
    gtmname << "force " << tgstr << ": total";
    gtm.stop(gtmname.str());

    logmsg.str("");
    logmsg << "Force assignment of " << required_per_petal << " "
        << tgstr << " targets";

    tm.stop();
    tm.report(logmsg.str().c_str());

    return;
}


void fba::Assignment::reassign_science_target(int32_t tstart, int32_t tstop,
    int32_t tile, int32_t fiber, int64_t target, bool balance_petals,
    std::map <int32_t, std::map <int32_t, bool> > const & done,
    int32_t & best_tile, int32_t & best_fiber) const {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    // Get the number of unused fibers on this current petal.
    int32_t petal = hw_->fiber_petal.at(fiber);
    int32_t passign = nassign_petal.at(TARGET_TYPE_SCIENCE).at(tile).at(petal);

    // Vector of available tile / fiber pairs which have a fiber that is a
    // science positioner.
    auto const & availtfall = favail_->data.at(target);
    std::vector < std::pair <int32_t, int32_t> > availtf;
    std::string pos_str("POS");
    for (auto const & av : availtfall) {
        if (pos_str.compare(hw_->fiber_device_type.at(av.second)) == 0) {
            availtf.push_back(av);
        }
    }

    best_tile = tile;
    best_fiber = fiber;
    int32_t best_passign = passign;

    if (extra_log) {
        logmsg.str("");
        logmsg << "reassign: tile " << tile << ", fiber "
            << fiber << ", target " << target
            << " considering for swap...";
        logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());

        logmsg.str("");
        logmsg << "reassign: tile " << tile << ", fiber "
            << fiber << ", target " << target << " considering tile indices "
            << tstart << " to " << tstop;
        logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
    }

    for (auto const & av : availtf) {
        // The available tile fiber pair
        int32_t av_tile = av.first;
        int32_t av_fiber = av.second;

        if (nassign_tile.at(TARGET_TYPE_SCIENCE).at(av_tile) == 0) {
            // This available tile / fiber is on a tile with
            // nothing assigned.  Skip it.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " available tile " << av_tile
                    << " has nothing assigned- skipping";
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            continue;
        }

        if ((done.count(av_tile) > 0)
            && (done.at(av_tile).count(av_fiber) > 0)
            && done.at(av_tile).at(av_fiber)) {
            // Already considered or swapped this available tile/fiber.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " avail T/F " << av_tile << "," << av_fiber
                    << " already considered or swapped";
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            continue;
        }

        if (fiber_target.at(av_tile).count(av_fiber) > 0) {
            // This available tile / fiber is already assigned.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " avail T/F " << av_tile << "," << av_fiber
                    << " already assigned";
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            continue;
        }

        if (target_fiber.at(target).count(av_tile) > 0) {
            // We have already assigned a fiber on this tile to this
            // target.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " already assigned on available tile " << av_tile;
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            continue;
        }

        int32_t av_petal = hw_->fiber_petal.at(av_fiber);
        int32_t av_tile_indx = tiles_->order.at(av_tile);

        if (av_tile_indx < tstart) {
            // The tile containing this alternate tile/fiber is before
            // the start of the tiles we are considering.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " available tile " << av_tile
                    << " at index " << av_tile_indx
                    << " is prior to tile start (" << tstart << ")";
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            continue;
        }

        // Projected target locations on the available tile.
        auto const & av_target_xy = tile_target_xy.at(av_tile);

        if ( ! ok_to_assign(hw_.get(), av_tile, av_fiber, target,
                            av_target_xy)) {
            // There must be a collision or some other problem.
            if (extra_log) {
                logmsg.str("");
                logmsg << "reassign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " avail T/F " << av_tile << "," << av_fiber
                    << " not OK to assign";
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            continue;
        }

        // Compare the number of assigned fibers on the petal of this
        // available tile/fiber to the current best one.
        int32_t av_passign =
            nassign_petal.at(TARGET_TYPE_SCIENCE).at(av_tile).at(av_petal);

        // At this point we know we have a new tile / fiber where we can place
        // this target.  If we are balancing the number of targets per petal,
        // check if this new tile / fiber is a better placement.  If we are not
        // balancing, then we just take this first available tile / fiber and
        // return.

        if (balance_petals) {
            if ((av_passign < hw_->nfiber_petal) &&
            (av_passign < best_passign)) {
                // There are some unassigned fibers on this available petal,
                // and the number of unassigned is greater than the original
                // tile/fiber.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "reassign: tile " << tile << ", fiber "
                        << fiber << ", target " << target
                        << " avail T/F " << av_tile << "," << av_fiber
                        << " has fewer assigned fibers in its petal ("
                        << av_passign << " vs "
                        << best_passign << "): new best";
                    logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
                }
                best_tile = av_tile;
                best_fiber = av_fiber;
                best_passign = av_passign;
            } else {
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "reassign: tile " << tile << ", fiber "
                        << fiber << ", target " << target
                        << " avail T/F " << av_tile << "," << av_fiber
                        << " has more assigned fibers in its petal ("
                        << av_passign << " vs "
                        << best_passign << "): skipping";
                    logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
                }
            }
        } else {
            best_tile = av_tile;
            best_fiber = av_fiber;
            break;
        }
    }

    return;
}


bool fba::Assignment::ok_to_assign (fba::Hardware const * hw, int32_t tile,
    int32_t fiber, int64_t target,
    std::map <int64_t, std::pair <double, double> > const & target_xy
    ) const {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    // Is the fiber stuck or broken?
    if (hw->state.at(fiber) != FIBER_STATE_OK) {
        if (extra_log) {
            logmsg.str("");
            logmsg << "ok_to_assign: tile " << tile << ", fiber "
                << fiber << " not OK";
            logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
        }
        return false;
    }

    if (fiber_target.count(tile) == 0) {
        // This is the first assignment on the tile, so definitely ok
        // to assign!
        if (extra_log) {
            logmsg.str("");
            logmsg << "ok_to_assign: tile " << tile
                << " first assignment (definitely OK)";
            logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
        }
        return true;
    }

    // Const reference to the assignment for this tile.
    auto const & ftile = fiber_target.at(tile);

    std::vector <int32_t> nbs;
    std::vector <int64_t> nbtarget;

    auto const & neighbors = hw->neighbors.at(fiber);

    // Check neighboring assignment.
    for (auto const & nb : neighbors) {
        if (ftile.count(nb) > 0) {
            // This neighbor has some assignment.
            int64_t nbtg = ftile.at(nb);
            if (nbtg == target) {
                // Target already assigned to a neighbor.
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "ok_to_assign: tile " << tile << ", fiber "
                        << fiber << ", target " << target
                        << " already assigned to neighbor fiber " << nb;
                    logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
                }
                return false;
            }
            nbs.push_back(nb);
            nbtarget.push_back(nbtg);
        }
    }

    // Would assigning this target produce a collision?  If a fiber is not
    // yet assigned, we assume that it is positioned at its central location
    // and so will not collide with anything.

    // Center position and target location for this fiber.
    fbg::dpair tcenter = hw->fiber_pos_xy_mm.at(fiber);
    fbg::dpair tpos = target_xy.at(target);

    size_t nnb = nbs.size();

    bool collide = false;

    // On average, the number of neighbors is 2-3.  Threading overhead seems
    // to negate the benefit here.
    // #pragma omp parallel for reduction(||:collide) schedule(dynamic) default(none) shared(hw, nbs, nbtarget, nnb, tcenter, tpos, target_xy)
    for (size_t b = 0; b < nnb; ++b) {
        int32_t const & nb = nbs[b];
        int64_t nbt = nbtarget[b];
        fbg::dpair ncenter;
        fbg::dpair npos;
        ncenter = hw->fiber_pos_xy_mm.at(nb);
        npos = target_xy.at(nbt);
        collide = hw->collide(tcenter, tpos, ncenter, npos);
        // Remove these lines if switching back to threading.
        if (collide) {
            if (extra_log) {
                logmsg.str("");
                logmsg << "ok_to_assign: tile " << tile << ", fiber "
                    << fiber << ", target " << target
                    << " would collide with target " << nbt;
                logger.debug_tfg(tile, fiber, target, logmsg.str().c_str());
            }
            return false;
        }
    }

    // Restore these lines if switching back to threading.
    // if (collide) {
    //     return false;
    // }

    // All good!
    return true;
}


int64_t fba::Assignment::find_best(fba::Hardware const * hw,
    fba::Targets * tgs, int32_t tile, int32_t fiber, uint8_t type,
    std::map <int64_t, std::pair <double, double> > const & target_xy,
    std::vector <int64_t> const & avail) const {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    bool extra_log = logger.extra_debug();

    int64_t best_id = -1;
    int32_t best_priority = -1;
    double best_subpriority = 0.0;
    int32_t best_obsremain = 0;

    if (avail.size() == 0) {
        // Nothing available
        return -1;
    }

    for (auto const & tgid : avail) {
        // This target
        auto const & tg = tgs->data.at(tgid);

        if ((type == TARGET_TYPE_SCIENCE) && (tg.obs_remain <= 0)) {
            // Skip science targets with no remaining observations.
            if (extra_log) {
                logmsg.str("");
                logmsg << "find_best: tile " << tile << ", fiber "
                    << fiber << ", target " << tgid
                    << " science target with no remaining obs";
                logger.debug_tfg(tile, fiber, tgid, logmsg.str().c_str());
            }
            continue;
        }

        // Is this target worth testing for assignment?
        bool test_target = false;

        if (type == TARGET_TYPE_SCIENCE) {
            // We are working with science targets.

            // FIXME:  This is where we could change the behavior of the way
            // that science targets are compared based on priority,
            // subpriority, and observations remaining.  The behavior below
            // replicates the selection of the original code.

            if (tg.priority > best_priority) {
                // Higher priority target class, always choose.
                test_target = true;
            } else if (tg.priority == best_priority) {
                // Same target class.  Use remaining obs and subpriority
                // to choose.
                if (tg.obs_remain > best_obsremain) {
                    test_target = true;
                } else if ( (tg.obs_remain == best_obsremain)
                            && (tg.subpriority > best_subpriority) ) {
                    test_target = true;
                }
            }
        } else {
            // We are working with either standards, skies or safe.  We base
            // our comparison on both priority and subpriority.  Normal
            // standards and skies will have priority == 0, so this reduces
            // to a comparison on subpriority.  Objects that are both a
            // standard and a science target will always "win", since their
            // priority is non-zero.  Change this behavior here if desired.
            double bestweight = static_cast <double> (best_priority)
                + best_subpriority;
            double newweight = static_cast <double> (tg.priority)
                + tg.subpriority;
            if (newweight > bestweight) {
                test_target = true;
            }
        }

        if (test_target) {
            if (extra_log) {
                logmsg.str("");
                logmsg << "find_best: tile " << tile << ", fiber "
                    << fiber << ", target " << tgid << ", type " << (int)tg.type
                    << " accept with priority = " << tg.priority
                    << ", subpriority = " << tg.subpriority
                    << ", obs_remain = " << tg.obs_remain;
                logger.debug_tfg(tile, fiber, tgid, logmsg.str().c_str());
            }
            if ((target_fiber.count(tgid) > 0) &&
                (target_fiber.at(tgid).count(tile) > 0)) {
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "find_best: tile " << tile << ", fiber "
                        << fiber << ", target " << tgid << ", type " << (int)tg.type
                        << " already assigned on current tile";
                    logger.debug_tfg(tile, fiber, tgid, logmsg.str().c_str());
                }
            } else if (ok_to_assign(hw, tile, fiber, tgid, target_xy)) {
                if (extra_log) {
                    logmsg.str("");
                    logmsg << "find_best: tile " << tile << ", fiber "
                        << fiber << ", target " << tgid << ", type "
                        << (int)tg.type << " SELECTED";
                    logger.debug_tfg(tile, fiber, tgid, logmsg.str().c_str());
                }
                best_id = tgid;
                best_priority = tg.priority;
                best_subpriority = tg.subpriority;
                best_obsremain = tg.obs_remain;
            }
        }
    }

    return best_id;
}


void fba::Assignment::assign_tilefiber(fba::Hardware const * hw,
    fba::Targets * tgs, int32_t tile, int32_t fiber, int64_t target,
    uint8_t type) {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    if (target < 0) {
        logmsg.str("");
        logmsg << "cannot assign negative target ID to tile "
            << tile << ", fiber " << fiber << ".  Did you mean to unassign?";
        logger.warning(logmsg.str().c_str());
        return;
    }

    int32_t petal = hw->fiber_petal.at(fiber);

    // if (fiber_target.count(tile) == 0) {
    //     fiber_target[tile].clear();
    // }
    auto & ftarg = fiber_target[tile];

    if (ftarg.count(fiber) > 0) {
        int64_t cur = ftarg.at(fiber);
        if (cur >= 0) {
            logmsg.str("");
            logmsg << "tile " << tile << ", fiber " << fiber
                << " already assigned to target " << cur
                << " cannot assign " << target;
            logger.warning(logmsg.str().c_str());
            return;
        }
    }

    // if (target_fiber.count(target) == 0) {
    //     target_fiber[target].clear();
    // }
    auto & tfiber = target_fiber[target];

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

    ftarg[fiber] = target;
    target_fiber[target][tile] = fiber;

    // Objects can be more than one type (e.g. standards and science).  When
    // incrementing the counts of object types per tile and petal, we want
    // to update the counts for valid types of this object.
    static const std::vector <uint8_t> target_types = {
        TARGET_TYPE_SCIENCE, TARGET_TYPE_STANDARD,
        TARGET_TYPE_SKY, TARGET_TYPE_SAFE};
    for (auto const & tt : target_types) {
        if (tgobj.is_type(tt)) {
            nassign_tile.at(tt).at(tile)++;
            nassign_petal.at(tt).at(tile).at(petal)++;
        }
    }
    tgobj.obs_remain--;

    return;
}


void fba::Assignment::unassign_tilefiber(fba::Hardware const * hw,
    fba::Targets * tgs, int32_t tile, int32_t fiber, uint8_t type) {

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    if (fiber_target.count(tile) == 0) {
        logmsg.str("");
        logmsg << "tile " << tile
            << " has no fibers assigned.  Ignoring unassign";
        logger.warning(logmsg.str().c_str());
        return;
    }
    auto & ftarg = fiber_target.at(tile);

    if (ftarg.count(fiber) == 0) {
        logmsg.str("");
        logmsg << "tile " << tile << ", fiber " << fiber
            << " already unassigned";
        logger.warning(logmsg.str().c_str());
        return;
    }

    int64_t target = ftarg.at(fiber);
    if (target < 0) {
        logmsg.str("");
        logmsg << "tile " << tile << ", fiber " << fiber
            << " already unassigned";
        logger.warning(logmsg.str().c_str());
        return;
    }

    auto & tgobj = tgs->data.at(target);

    if ( ! tgobj.is_type(type)) {
        logmsg.str("");
        logmsg << "current target " << target << " not of type "
            << (int)type << " requested in unassign of tile " << tile
            << ", fiber " << fiber;
        logger.error(logmsg.str().c_str());
        throw std::runtime_error(logmsg.str().c_str());
    }

    int32_t petal = hw->fiber_petal.at(fiber);

    // Objects can be more than one type (e.g. standards and science).  When
    // incrementing the counts of object types per tile and petal, we want
    // to update the counts for valid types of this object.
    static const std::vector <uint8_t> target_types = {
        TARGET_TYPE_SCIENCE, TARGET_TYPE_STANDARD,
        TARGET_TYPE_SKY, TARGET_TYPE_SAFE};
    for (auto const & tt : target_types) {
        if (tgobj.is_type(tt)) {
            nassign_tile.at(tt).at(tile)--;
            nassign_petal.at(tt).at(tile).at(petal)--;
        }
    }
    tgobj.obs_remain++;

    target_fiber[target].erase(tile);
    ftarg.erase(fiber);

    return;
}


void fba::Assignment::targets_to_project(
    fba::Targets const * tgs,
    std::map <int32_t, std::vector <int64_t> > const & tgsavail,
    std::vector <int32_t> const & fibers,
    std::vector <int64_t> & tgids,
    std::vector <double> & tgra,
    std::vector <double> & tgdec) const {
    // This function computes the target IDs that need to be projected
    // for a given set of fibers on a tile.

    std::set <int64_t> seen;
    tgids.clear();
    tgra.clear();
    tgdec.clear();

    for (auto const & fid : fibers) {
        // Project this fiber's targets
        if (tgsavail.count(fid) > 0) {
            // The available targets for this fiber.
            auto const & avail = tgsavail.at(fid);
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
        int32_t tile_id, double tile_ra, double tile_dec,
        std::map <int64_t, std::pair <double, double> > & target_xy) const {
    // This function computes the projection of all available targets
    // into focalplane coordinates for one tile.  We do this once and then
    // use these positions multiple times.

    // Clear the output
    target_xy.clear();

    if (tgsavail->data.count(tile_id) == 0) {
        // This tile has no fibers with available targets.
        return;
    }

    // Reference to the available targets for this tile.
    auto const & tfavail = tgsavail->data.at(tile_id);

    if (tfavail.size() == 0) {
        // There are no fibers on this tile with available targets.
        return;
    }

    // List of fibers we are projecting- anything with targets available.
    std::vector <int32_t> fids;
    for (auto const & f : hw->fiber_id) {
        if ((tfavail.count(f) != 0) && (tfavail.at(f).size() > 0)) {
            fids.push_back(f);
        }
    }

    // Vectors of all targets to compute
    std::vector <int64_t> tgids;
    std::vector <double> tgra;
    std::vector <double> tgdec;

    targets_to_project(tgs, tfavail, fids, tgids, tgra, tgdec);

    // Now thread over the targets to compute

    std::vector <std::pair <double, double> > xy;
    hw->radec2xy_multi(tile_ra, tile_dec, tgra, tgdec, xy);

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


fba::FibersAvailable::pshr fba::Assignment::fibers_avail() const {
    return favail_;
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
