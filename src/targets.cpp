// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <cmath>
#include <cstdint>
#include <limits>

#include <utils.h>
#include <tiles.h>
#include <targets.h>
#include <assert.h>

namespace fba = fiberassign;

namespace fbg = fiberassign::geom;


std::string fba::target_string(uint8_t type) {
    std::string ret;
    if (type == TARGET_TYPE_SCIENCE) {
        ret = "science";
    } else if (type == TARGET_TYPE_SKY) {
        ret = "sky";
    } else if (type == TARGET_TYPE_STANDARD) {
        ret = "standard";
    } else if (type == TARGET_TYPE_SAFE) {
        ret = "safe";
    } else if (type == TARGET_TYPE_SUPPSKY) {
        ret = "suppsky";
    } else {
        ret = "NA";
        throw std::runtime_error("Unknown target type");
    }
    return ret;
}


fba::Target::Target() {
    id = -1;
    obsremain = 0;
    priority = 0;
    subpriority = 0.0;
    obscond = 0;
    type = 0;
}


fba::Target::Target(int64_t tid,
                    int32_t tobsremain,
                    int32_t tpriority,
                    double tsubpriority,
                    int32_t tobscond,
                    uint8_t ttype) {
    id = tid;
    obsremain = tobsremain;
    priority = tpriority;
    subpriority = tsubpriority;
    obscond = tobscond;
    type = ttype;
}


bool fba::Target::is_science() const {
    return ((type & TARGET_TYPE_SCIENCE) != 0);
}


bool fba::Target::is_standard() const {
    return ((type & TARGET_TYPE_STANDARD) != 0);
}


bool fba::Target::is_sky() const {
    return ((type & TARGET_TYPE_SKY) != 0);
}


bool fba::Target::is_suppsky() const {
    return ((type & TARGET_TYPE_SUPPSKY) != 0);
}


bool fba::Target::is_safe() const {
    return ((type & TARGET_TYPE_SAFE) != 0);
}


bool fba::Target::is_type(uint8_t t) const {
    return ((type & t) != 0);
}


double fba::Target::total_priority() const {
    // This is where we could control the depth-first vs. breadth-first
    // behavior.  Default is breadth-first:
    return (double)(priority * 100 + obsremain) + subpriority;
    //
    // Instead, we probably want to use some bit in the target bits to
    // select whether to prioritize depth vs breadth first.
    // For example:
    // if (bits & DEPTH_FIRST_MASK) {
    //     return (double)(priority * 100 + (100 - obsremain)) + subpriority;
    // } else {
    //     return (double)(priority * 100 + obsremain) + subpriority;
    // }
}


fba::Targets::Targets() {
    science_classes.clear();
    survey = "";
}


void fba::Targets::append(std::string const & tsurvey,
                          std::vector <int64_t> const & id,
                          std::vector <int32_t> const & obsremain,
                          std::vector <int32_t> const & priority,
                          std::vector <double> const & subpriority,
                          std::vector <int32_t> const & obscond,
                          std::vector <uint8_t> const & type) {
    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    if (survey.compare("") == 0) {
        survey = tsurvey;
    } else if (survey.compare(tsurvey) != 0) {
        logmsg.str("");
        logmsg << "Targets object has survey type \"" << survey
            << "\", cannot append data from survey type \""
            << tsurvey << "\"";
        logger.error(logmsg.str().c_str());
        throw std::runtime_error(logmsg.str().c_str());
    }

    for (size_t t = 0; t < id.size(); ++t) {
        if (type[t] == 0) {
            // This target is not one of the recognized categories (science,
            // standard, sky, suppsky, or safe).  Skip it.
            logmsg.str("");
            logmsg << "Survey " << survey
                << " target ID " << id[t]
                << " type not identified- SKIPPING";
            logger.debug(logmsg.str().c_str());
            continue;
        }
        if (data.count(id[t]) > 0) {
            // This target already exists.  This is an error.
            logmsg.str("");
            auto const & tg = data.at(id[t]);
            logmsg << "Target ID " << id[t]
                << " already exists with properties: (" << tg.priority << ","
                << tg.subpriority << ") " << tg.obsremain << ", "
                << tg.obscond << ", " << (int)(tg.type);
            logger.error(logmsg.str().c_str());
            logmsg.str("");
            logmsg << "  New target properties: (" << priority[t] << ","
                << subpriority[t] << ") " << obsremain[t] << ", "
                << obscond[t];
            logger.error(logmsg.str().c_str());
            logmsg.str("");
            logmsg << "  Duplicate target IDs are not permitted.";
            logger.error(logmsg.str().c_str());
            throw std::runtime_error(logmsg.str().c_str());
        } else {
            data[id[t]] = Target(id[t],
                                 obsremain[t], priority[t], subpriority[t],
                                 obscond[t], type[t]);
        }
        if ((priority[t] > 0) && ((type[t] & TARGET_TYPE_SCIENCE) != 0)) {
            // Only consider science targets in the list of target classes.
            if (science_classes.count(priority[t]) == 0) {
                science_classes.insert(priority[t]);
            }
        }
    }
    // for (auto const & tc : science_classes) {
    //     logmsg.str("");
    //     logmsg << "Targets:  have science target class " << tc;
    //     logger.debug(logmsg.str().c_str());
    // }
}


fba::TargetsAvailable::TargetsAvailable(Hardware::pshr hw,
                                        Tiles::pshr tiles,
                                        std::map<int64_t, std::vector<int64_t> > tile_targetids,
                                        std::map<int64_t, std::vector<double> > tile_x,
                                        std::map<int64_t, std::vector<double> > tile_y) {
    Timer tm;
    tm.start();

    fba::Logger & logger = fba::Logger::get();

    data.clear();

    tiles_ = tiles;
    hw_ = hw;

    size_t ntile = tiles_->id.size();
    size_t nloc = hw_->nloc;

    // patrol buffer
    double patrol_buffer = hw_->patrol_buffer_mm;

    std::vector <int32_t> loc(nloc);
    std::vector <double> loc_center_x(nloc);
    std::vector <double> loc_center_y(nloc);
    std::vector <double> loc_patrol(nloc);

    for (size_t j = 0; j < nloc; ++j) {
        loc[j] = hw_->locations[j];
        loc_center_x[j] = hw_->loc_pos_curved_mm[loc[j]].first;
        loc_center_y[j] = hw_->loc_pos_curved_mm[loc[j]].second;
        loc_patrol[j] = hw_->loc_theta_arm[loc[j]] + hw_->loc_phi_arm[loc[j]]
            - patrol_buffer;
    }

    // shared_ptr reference counting is not threadsafe.  Here we extract
    // a copy of the "raw" pointers needed inside the parallel region.

    Tiles * ptiles = tiles.get();
    Hardware * phw = hw_.get();

    #pragma omp parallel default(shared)
    {
        // We re-use these thread-local vectors to reduce the number
        // of times we are realloc'ing memory for every fiber.
        std::vector <int64_t> nearby;
        std::vector <KdTreePoint> nearby_tree_points;
        std::vector <KdTreePoint> nearby_data;

        std::ostringstream logmsg;

        // Thread local properties.
        int32_t tid;

        // Thread local data for available targets- reduced at the end.
        std::map <int32_t, std::map <int32_t,
                                     std::vector <int64_t> > > thread_data;
        std::map <int32_t, std::map <int32_t,
                                     std::vector <
                                         std::pair<double, double> > > > thread_data_xy;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < ntile; ++i) {
            tid = ptiles->id[i];
            thread_data[tid].clear();
            thread_data_xy[tid].clear();

            if (tile_targetids[tid].size() == 0) {
                // No targets for this tile.
                continue;
            }

            assert(tile_targetids[tid].size() == tile_x[tid].size());
            assert(tile_targetids[tid].size() == tile_y[tid].size());

            nearby_tree_points.clear();
            KdTreePoint cur;
            auto vx = tile_x[tid].begin();
            auto vy = tile_y[tid].begin();
            for (auto vid = tile_targetids[tid].begin();
                 vid != tile_targetids[tid].end(); vid++, vx++, vy++) {
                cur.id = *vid;
                cur.pos[0] = *vx;
                cur.pos[1] = *vy;
                nearby_tree_points.push_back(cur);
            }

            KDtree <KdTreePoint> nearby_tree(nearby_tree_points, 2);

            size_t locs_with_targets = 0;

            double loc_pos[2];

            for (size_t j = 0; j < nloc; ++j) {
                thread_data[tid][loc[j]].resize(0);
                thread_data_xy[tid][loc[j]].resize(0);
                loc_pos[0] = loc_center_x[j];
                loc_pos[1] = loc_center_y[j];

                // Lookup targets near this location in focalplane
                // coordinates.
                nearby_data = nearby_tree.near_with_data(loc_pos, 0.0, loc_patrol[j]);
                if (nearby_data.size() == 0) {
                    // No targets for this location
                    continue;
                }

                // The kdtree gets us the targets that are close to
                // our region of interest around each positioner.  Now
                // go through these targets and check whether the
                // positioner can move to each one.  We DO NOT sort
                // targets by priority, since the total priority will
                // change as observations are made.  Instead, this
                // sorting is done for each tile during assignment.

                for (auto const & tnear : nearby_data) {
                    fbg::dpair obj_xy;
                    obj_xy.first  = tnear.pos[0];
                    obj_xy.second = tnear.pos[1];
                    bool fail = phw->position_xy_bad(loc[j], obj_xy);
                    if (fail) {
                        if (logger.extra_debug()) {
                            logmsg.str("");
                            logmsg << std::setprecision(2) << std::fixed;
                            logmsg << "targets avail:  tile " << tid
                                << ", loc " << loc[j] << ", kdtree target "
                                << tnear.id << " not physically reachable by positioner";
                            logger.debug_tfg(tid, loc[j], tnear.id,
                                             logmsg.str().c_str());
                        }
                    } else {
                        thread_data[tid][loc[j]].push_back(tnear.id);
                        thread_data_xy[tid][loc[j]].push_back(std::make_pair(tnear.pos[0], tnear.pos[1]));
                    }
                }

                if (thread_data[tid][loc[j]].size() > 0) {
                    locs_with_targets++;
                }
            }
        }

        // Now reduce across threads.
        #pragma omp critical
        {
            for (auto const & it : thread_data) {
                int32_t ttile = it.first;

                if (data.count(ttile) > 0) {
                    throw std::runtime_error("Available target data already exists for tile");
                }

                data[ttile].clear();

                auto const & fmap = it.second;

                for (size_t f = 0; f < nloc; ++f) {
                    int32_t lid = loc[f];
                    if (fmap.count(lid) == 0) {
                        // This location had no targets.
                        continue;
                    }
                    data[ttile][lid] = fmap.at(lid);
                }
            }
            thread_data.clear();
            for (auto const & it : thread_data_xy) {
                int32_t ttile = it.first;
                if (data_xy.count(ttile) > 0) {
                    throw std::runtime_error("Available target data already exists for tile");
                }
                data_xy[ttile].clear();
                auto const & fmap = it.second;
                for (size_t f = 0; f < nloc; ++f) {
                    int32_t lid = loc[f];
                    if (fmap.count(lid) == 0) {
                        // This location had no targets.
                        continue;
                    }
                    data_xy[ttile][lid] = fmap.at(lid);
                }
            }
            thread_data_xy.clear();
        }
    }

    for (size_t i = 0; i < ntile; ++i) {
        int32_t tid = ptiles->id[i];
        if (data.count(tid) == 0) {
            continue;
        }
        size_t total_avail = 0;
        for (size_t j = 0; j < nloc; ++j) {
            if (data.at(tid).count(loc[j]) > 0) {
                total_avail += data.at(tid).at(loc[j]).size();
            }
        }
        std::ostringstream msg;
        msg.str("");
        msg << "targets avail:  tile " << tid
            << ", " << total_avail << " total available targets";
        logger.debug(msg.str().c_str());
    }

    tm.stop();
    tm.report("Computing targets available to all tile / locations");
}

fba::Hardware::pshr fba::TargetsAvailable::hardware() const {
    return hw_;
}

fba::Tiles::pshr fba::TargetsAvailable::tiles() const {
    return tiles_;
}


std::map <int32_t, std::vector <int64_t> > fba::TargetsAvailable::tile_data(int32_t tile) const {
    if (data.count(tile) == 0) {
        return std::map <int32_t, std::vector <int64_t> >();
    } else {
        return data.at(tile);
    }
}


fba::LocationsAvailable::LocationsAvailable(fba::TargetsAvailable::pshr tgsavail) {
    fba::Timer tm;
    tm.start();

    fba::Logger & logger = fba::Logger::get();

    data.clear();

    // In order to play well with OpenMP for loops, construct a simple vector of
    // std::map keys for the available objects.
    std::vector <int32_t> tfkeys;
    for (auto const & tfit : tgsavail->data) {
        tfkeys.push_back(tfit.first);
    }

    size_t ntile = tfkeys.size();

    // shared_ptr reference counting is not threadsafe.  Here we extract
    // a copy of the "raw" pointers needed inside the parallel region.

    auto * ptgsavail = tgsavail.get();

    #pragma omp parallel default(none) shared(logger, ptgsavail, ntile, tfkeys)
    {
        // Our thread-local data, to be reduced at the end.
        std::map < int64_t,
            std::vector < std::pair <int32_t, int32_t> > > thread_data;
        std::ostringstream logmsg;

        auto const & avail = ptgsavail->data;

        #pragma omp for schedule(dynamic)
        for (size_t tindx = 0; tindx < ntile; ++tindx) {
            int32_t tid = tfkeys[tindx];
            if (avail.count(tid) == 0) {
                continue;
            }
            auto tavail = avail.at(tid);
            for (auto const & ltg : tavail) {
                auto loc = ltg.first;
                for (auto const & tg : ltg.second) {
                    if (thread_data.count(tg) == 0) {
                        thread_data[tg].resize(0);
                    }
                    if (logger.extra_debug()) {
                        logmsg.str("");
                        logmsg << "target " << tg
                            << " has available tile / location "
                            << tid << ", " << loc;
                        logger.debug_tfg(tid, loc, tg, logmsg.str().c_str());
                    }
                    thread_data[tg].push_back(std::make_pair(tid, loc));
                }
            }
        }

        // Now reduce across threads.
        #pragma omp critical
        {
            for (auto const & it : thread_data) {
                int64_t tg = it.first;
                if (data.count(tg) == 0) {
                    data[tg].resize(0);
                }
                for (auto const & ft : it.second) {
                    data[tg].push_back(ft);
                }
            }
            thread_data.clear();
        }
    }

    // Sort the available tile / locations by tile order, since the original
    // order will depend on thread concurrency.

    auto tls = tgsavail->tiles();
    auto torder = tls->order;
    auto tids = tls->id;

    fba::tile_loc_compare tlcomp;

    std::vector <tile_loc> temptl;

    for (auto & tgav : data) {
        auto & av = tgav.second;
        temptl.clear();
        for (auto const & tf : av) {
            temptl.push_back(
                std::make_pair(torder.at(tf.first), tf.second)
            );
        }
        std::stable_sort(temptl.begin(), temptl.end(), tlcomp);
        av.clear();
        for (auto const & tf : temptl) {
            av.push_back(
                std::make_pair(tids.at(tf.first), tf.second)
            );
        }
    }

    tm.stop();
    tm.report("Computing tile / locations available to all objects");
}


std::vector <std::pair <int32_t, int32_t> >
    fba::LocationsAvailable::target_data(int64_t target) const {
    if (data.count(target) == 0) {
        return std::vector <std::pair <int32_t, int32_t> >();
    } else {
        return data.at(target);
    }
}
