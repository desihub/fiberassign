// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <cmath>
#include <cstdint>
#include <limits>

#include <utils.h>
#include <tiles.h>
#include <targets.h>


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
    } else {
        ret = "NA";
        throw std::runtime_error("Unknown target type");
    }
    return ret;
}


fba::Target::Target() {
    id = -1;
    ra = 0.0;
    dec = 0.0;
    bits = 0;
    obsremain = 0;
    priority = 0;
    subpriority = 0.0;
    obscond = 0;
    type = 0;
}


fba::Target::Target(int64_t tid,
                    double tra,
                    double tdec,
                    int64_t tbits,
                    int32_t tobsremain,
                    int32_t tpriority,
                    double tsubpriority,
                    int32_t tobscond,
                    uint8_t ttype) {
    id = tid;
    ra = tra;
    dec = tdec;
    bits = tbits;
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


bool fba::Target::is_safe() const {
    return ((type & TARGET_TYPE_SAFE) != 0);
}


bool fba::Target::is_type(uint8_t t) const {
    return ((type & t) != 0);
}


fba::Targets::Targets() {
    science_classes.clear();
    survey = "";
}


void fba::Targets::append(std::string const & tsurvey,
                          std::vector <int64_t> const & id,
                          std::vector <double> const & ra,
                          std::vector <double> const & dec,
                          std::vector <int64_t> const & targetbits,
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
            // standard, sky, or safe).  Skip it.
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
                << " already exists with properties: ("
                << tg.ra << "," << tg.dec << ") (" << tg.priority << ","
                << tg.subpriority << ") " << tg.obsremain << ", "
                << tg.obscond << ", " << (int)(tg.type);
            logger.error(logmsg.str().c_str());
            logmsg.str("");
            logmsg << "  New target properties: ("
                << ra[t] << "," << dec[t] << ") (" << priority[t] << ","
                << subpriority[t] << ") " << obsremain[t] << ", "
                << obscond[t];
            logger.error(logmsg.str().c_str());
            logmsg.str("");
            logmsg << "  Duplicate target IDs are not permitted.";
            logger.error(logmsg.str().c_str());
            throw std::runtime_error(logmsg.str().c_str());
        } else {
            data[id[t]] = Target(id[t], ra[t], dec[t], targetbits[t],
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


fba::TargetTree::TargetTree(Targets::pshr objs, double min_tree_size) {
    Timer tm;
    tm.start();

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    mintreesz_ = min_tree_size;

    treelist_.resize(0);

    TreePoint tp;
    double theta;
    double phi;
    double deg2rad = M_PI / 180.0;
    double stheta;

    for (auto const & obj : objs->data) {
        tp.id = obj.first;
        theta = (90.0 - obj.second.dec) * deg2rad;
        phi = (obj.second.ra) * deg2rad;
        stheta = ::sin(theta);
        tp.nhat[0] = ::cos(phi) * stheta;
        tp.nhat[1] = ::sin(phi) * stheta;
        tp.nhat[2] = ::cos(theta);
        if (logger.extra_debug()) {
            logmsg.str("");
            logmsg << "add target ID " << obj.first << " to tree at "
                << "RA = " << obj.second.ra << ", DEC = " << obj.second.dec;
            logger.debug_tfg(-1, -1, obj.first, logmsg.str().c_str());
        }
        treelist_.push_back(tp);
    }

    tree_ = std::make_unique < htmTree <TreePoint> > (treelist_, mintreesz_);
    tree_->stats();
    tm.stop();
    tm.report("Building target tree");
}


void fba::TargetTree::near(double ra_deg, double dec_deg, double radius_rad,
    std::vector <int64_t> & result) const {

    double deg2rad = M_PI / 180.0;
    double theta = (90.0 - dec_deg) * deg2rad;
    double phi = ra_deg * deg2rad;
    double stheta = ::sin(theta);
    double nhat[3];
    nhat[0] = ::cos(phi) * stheta;
    nhat[1] = ::sin(phi) * stheta;
    nhat[2] = ::cos(theta);

    std::vector <int> temp = tree_->near (treelist_, nhat, radius_rad);
    result.resize(temp.size());
    for (size_t i = 0; i < temp.size(); ++i) {
        result[i] = treelist_[temp[i]].id;
    }
    return;
}


fba::TargetsAvailable::TargetsAvailable(Hardware::pshr hw, Targets::pshr objs,
                                        Tiles::pshr tiles,
                                        TargetTree::pshr tree) {
    Timer tm;
    tm.start();

    fba::Logger & logger = fba::Logger::get();

    data.clear();

    double deg2rad = M_PI / 180.0;

    tiles_ = tiles;
    hw_ = hw;

    size_t ntile = tiles_->id.size();
    size_t nloc = hw_->nloc;

    // Radius of the tile.
    double tile_radius = hw_->focalplane_radius_deg * deg2rad;

    std::vector <int32_t> loc(nloc);
    std::vector <double> loc_center_x(nloc);
    std::vector <double> loc_center_y(nloc);
    std::vector <double> loc_patrol(nloc);

    for (size_t j = 0; j < nloc; ++j) {
        loc[j] = hw_->locations[j];
        double cx = hw_->loc_pos_xy_mm[loc[j]].first;
        double cy = hw_->loc_pos_xy_mm[loc[j]].second;
        loc_center_x[j] = cx;
        loc_center_y[j] = cy;

        // FIXME:  This line only exists to achieve consistency with the
        // legacy results when using a focalplane model based on a fake
        // fiberpos file.  Remove this after this version is tagged.
        loc_patrol[j] = 5.8;
        // loc_patrol[j] = hw_->loc_theta_arm[loc[j]] + hw_->loc_phi_arm[loc[j]];
    }

    // shared_ptr reference counting is not threadsafe.  Here we extract
    // a copy of the "raw" pointers needed inside the parallel region.

    Targets * pobjs = objs.get();
    Tiles * ptiles = tiles.get();
    TargetTree * ptree = tree.get();
    Hardware * phw = hw_.get();

    #pragma omp parallel default(none) shared(logger, pobjs, ptiles, ptree, phw, ntile, tile_radius, nloc, loc, loc_center_x, loc_center_y, loc_patrol)
    {
        // We re-use these thread-local vectors to reduce the number
        // of times we are realloc'ing memory for every fiber.
        std::vector <int64_t> nearby;
        std::vector <KdTreePoint> nearby_tree_points;
        std::vector <int64_t> nearby_loc;
        std::vector <int64_t> result;
        std::vector <weight_index> result_weight;
        std::ostringstream logmsg;

        // Thread local properties.
        int32_t tid;
        double tra;
        double tdec;
        int32_t tobs;
        std::pair <double, double> target_xy;
        weight_compare result_comp;

        // Thread local data for available targets- reduced at the end.
        std::map <int32_t, std::map <int32_t,
                  std::vector <int64_t> > > thread_data;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < ntile; ++i) {
            tid = ptiles->id[i];
            tra = ptiles->ra[i];
            tdec = ptiles->dec[i];
            tobs = ptiles->obscond[i];
            thread_data[tid].clear();

            // Get all targets in range of this tile
            ptree->near(tra, tdec, tile_radius, nearby);

            if (nearby.size() == 0) {
                // No targets for this tile
                continue;
            }

            // Project objects onto the focal plane

            nearby_tree_points.clear();
            KdTreePoint cur;
            for (auto const & tnear : nearby) {
                auto const & obj = pobjs->data[tnear];
                if ((obj.obscond & tobs) != 0) {
                    // Only use targets with correct obs conditions
                    target_xy = phw->radec2xy(tra, tdec, obj.ra, obj.dec);
                    cur.id = obj.id;
                    cur.pos[0] = target_xy.first;
                    cur.pos[1] = target_xy.second;
                    nearby_tree_points.push_back(cur);
                } else if (logger.extra_debug()) {
                    logmsg.str("");
                    logmsg << std::setprecision(2) << std::fixed;
                    logmsg << "targets avail:  tile " << tid
                        << ", skip TARGET ID " << obj.id
                        << " with incompatible OBSCONDITIONS ("
                        << obj.obscond << ") for tile (" << tobs << ")";
                    logger.debug_tfg(tid, -1, obj.id,
                                     logmsg.str().c_str());
                }
            }
            KDtree <KdTreePoint> nearby_tree(nearby_tree_points, 2);

            size_t locs_with_targets = 0;

            double loc_pos[2];

            for (size_t j = 0; j < nloc; ++j) {
                thread_data[tid][loc[j]].resize(0);
                loc_pos[0] = loc_center_x[j];
                loc_pos[1] = loc_center_y[j];
                auto loc_xy = std::make_pair(loc_center_x[j],
                                             loc_center_y[j]);

                // Lookup targets near this location in focalplane
                // coordinates.
                nearby_loc = nearby_tree.near(loc_pos, 0.0, loc_patrol[j]);

                if (nearby_loc.size() == 0) {
                    // No targets for this location
                    continue;
                }

                // The kdtree gets us the targets that are close to our
                // region of interest.  Now go through these targets and
                // compute the precise distance from the loc center.
                // We sort the available targets by priority now,
                // to avoid doing it multiple times later.
                result.clear();
                result_weight.clear();

                size_t tindx = 0;
                for (auto const & tnear : nearby_loc) {
                    auto & obj = pobjs->data[tnear];
                    auto obj_xy = phw->radec2xy(tra, tdec, obj.ra,
                                                obj.dec);
                    double dist = geom::sq(loc_xy, obj_xy);
                    double pdist = geom::sq(loc_patrol[j]);
                    if (dist > pdist) {
                        // outside the patrol radius
                        if (logger.extra_debug()) {
                            logmsg.str("");
                            logmsg << std::setprecision(2) << std::fixed;
                            logmsg << "targets avail:  tile " << tid
                                << ", loc " << loc[j] << ", kdtree target "
                                << obj.id << " outside patrol radius ("
                                << dist << " > " << pdist << ")";
                            logger.debug_tfg(tid, loc[j], obj.id,
                                             logmsg.str().c_str());
                        }
                        continue;
                    }
                    double totpriority =
                        static_cast <double> (obj.priority)
                        + obj.subpriority;
                    result.push_back(tnear);
                    result_weight.push_back(
                        std::make_pair(totpriority, tindx));
                    tindx++;
                }

                std::stable_sort(result_weight.begin(),
                                 result_weight.end(), result_comp);

                for (auto const & wt : result_weight) {
                    if (logger.extra_debug()) {
                        logmsg.str("");
                        logmsg << std::setprecision(2) << std::fixed;
                        logmsg << "targets avail:  tile " << tid
                            << ", location " << loc[j]
                            << " append ID "
                            << result[wt.second] << " (type="
                            << (int)(pobjs->data[result[wt.second]].type) << ")"
                            << ", total priority " << wt.first;
                        logger.debug_tfg(tid, loc[j],
                                         result[wt.second],
                                         logmsg.str().c_str());
                    }
                    thread_data[tid][loc[j]].push_back(
                        result[wt.second]);
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

    //std::cout << "LocationsAvailable:  tgsavail has " << avail.size() << " tiles" << std::endl;

    // In order to play well with OpenMP for loops, construct a simple vector of
    // std::map keys for the available objects.
    std::vector <int32_t> tfkeys;
    for (auto const & tfit : tgsavail->data) {
        tfkeys.push_back(tfit.first);
    }

    size_t ntile = tfkeys.size();

    // This structure is only used in order to reproduce the previous
    // behavior of looping in fiber ID order.
    auto loc = tgsavail->hardware()->locations;
    auto loc_to_fiber = tgsavail->hardware()->loc_fiber;
    std::map <int32_t, int32_t> fiber_to_loc;
    std::vector <fba::fiber_loc> fiber_and_loc;
    for (auto const & lid : loc) {
        fiber_to_loc[loc_to_fiber[lid]] = lid;
        fiber_and_loc.push_back(std::make_pair(loc_to_fiber[lid], lid));
    }
    fba::fiber_loc_compare flcomp;
    std::stable_sort(fiber_and_loc.begin(), fiber_and_loc.end(), flcomp);

    // shared_ptr reference counting is not threadsafe.  Here we extract
    // a copy of the "raw" pointers needed inside the parallel region.

    auto * ptgsavail = tgsavail.get();

    #pragma omp parallel default(none) shared(logger, ptgsavail, ntile, tfkeys, fiber_and_loc)
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
            for (auto const & flid : fiber_and_loc) {
                auto fid = flid.first;
                auto lid = flid.second;
                if (tavail.count(lid) == 0) {
                    continue;
                }
                auto lavail = tavail.at(lid);
                for (auto const & tg : lavail) {
                    if (thread_data.count(tg) == 0) {
                        thread_data[tg].resize(0);
                    }
                    if (logger.extra_debug()) {
                        logmsg.str("");
                        logmsg << "target " << tg
                            << " has available tile / location "
                            << tid << ", " << lid
                            << " (fiber " << fid << ")";
                        logger.debug_tfg(tid, lid, tg, logmsg.str().c_str());
                    }
                    thread_data[tg].push_back(std::make_pair(tid, lid));
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

    std::vector < std::pair <int32_t, int32_t> > tdx;

    for (auto & tgav : data) {
        // int64_t tgid = tgav.first;
        auto & av = tgav.second;
        tdx.clear();
        for (auto const & tf : av) {
            tdx.push_back(std::make_pair(torder.at(tf.first),
                          loc_to_fiber.at(tf.second)));
        }
        std::stable_sort(tdx.begin(), tdx.end(), tlcomp);
        av.clear();
        for (auto const & tf : tdx) {
            av.push_back(std::make_pair(tids[tf.first],
                         fiber_to_loc.at(tf.second)));
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
