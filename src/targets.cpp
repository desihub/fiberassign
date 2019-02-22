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
    desi_target = 0;
    bgs_target = 0;
    mws_target = 0;
    obs_remain = 0;
    priority = 0;
    subpriority = 0.0;
    obscond = 0;
    type = 0;
}


fba::Target::Target(int64_t tid,
                    double tra,
                    double tdec,
                    int64_t tdesi_target,
                    int64_t tbgs_target,
                    int64_t tmws_target,
                    int32_t tobs_remain,
                    int32_t tpriority,
                    double tsubpriority,
                    int32_t tobscond,
                    uint8_t ttype) {
    id = tid;
    ra = tra;
    dec = tdec;
    desi_target = tdesi_target;
    bgs_target = tbgs_target;
    mws_target = tmws_target;
    obs_remain = tobs_remain;
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
}


void fba::Targets::append(std::vector <int64_t> const & id,
                          std::vector <double> const & ra,
                          std::vector <double> const & dec,
                          std::vector <int64_t> const & desi_target,
                          std::vector <int64_t> const & bgs_target,
                          std::vector <int64_t> const & mws_target,
                          std::vector <int32_t> const & obs_remain,
                          std::vector <int32_t> const & priority,
                          std::vector <double> const & subpriority,
                          std::vector <int32_t> const & obscond,
                          std::vector <uint8_t> const & type) {
    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;
    for (size_t t = 0; t < id.size(); ++t) {
        if (type[t] == 0) {
            // This target is not one of the recognized categories (science,
            // standard, sky, or safe).  Skip it.
            continue;
        }
        if (data.count(id[t]) > 0) {
            // This target already exists.  Check that its location and
            // properties match *except* for its type.  The object might
            // be specified twice with different bits for the type, so we
            // OR these together.
            auto & tg = data.at(id[t]);
            bool dup = true;
            if (::fabs((tg.ra - ra[t])/tg.ra)
                > std::numeric_limits<float>::epsilon()) {
                dup = false;
            }
            if (::fabs((tg.dec - dec[t])/tg.dec)
                > std::numeric_limits<float>::epsilon()) {
                dup = false;
            }
            if (dup) {
                // bitwise or the type
                tg.type |= type[t];
                // bitwise or the obs cond
                tg.obscond |= obscond[t];
                // bitwise or the target masks
                tg.desi_target |= desi_target[t];
                tg.bgs_target |= bgs_target[t];
                tg.mws_target |= mws_target[t];
                // choose the larger of the obs remaining
                if (obs_remain[t] > tg.obs_remain) {
                    tg.obs_remain = obs_remain[t];
                }
                // choose the larger of the priority and subpriority
                if (priority[t] > tg.priority) {
                    tg.priority = priority[t];
                }
                if (subpriority[t] > tg.subpriority) {
                    tg.subpriority = subpriority[t];
                }
            } else {
                // this is a conflicting object
                logmsg.str("");
                logmsg << "Target ID " << id[t]
                    << " already exists with properties: ("
                    << tg.ra << "," << tg.dec << ") (" << tg.priority << "," << tg.subpriority << ") " << tg.obs_remain << ", " << tg.obscond << ", " << (int)(tg.type);
                logger.warning(logmsg.str().c_str());
                logmsg.str("");
                logmsg << "  New target properties: ("
                    << ra[t] << "," << dec[t] << ") (" << priority[t] << "," << subpriority[t] << ") " << obs_remain[t] << ", " << obscond[t];
                logger.warning(logmsg.str().c_str());
                logmsg.str("");
                logmsg << "  Ignoring new target info";
                logger.warning(logmsg.str().c_str());
            }
        } else {
            data[id[t]] = Target(id[t], ra[t], dec[t], desi_target[t],
                                 bgs_target[t], mws_target[t], obs_remain[t],
                                 priority[t], subpriority[t], obscond[t],
                                 type[t]);
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
    size_t nfiber = hw_->nfiber;

    // Radius of the tile.
    double tile_radius = hw_->focalplane_radius_deg * deg2rad;

    // Patrol radius in mm on focalplane.
    double fiber_patrol_mm = hw_->patrol_mm;

    std::vector <int32_t> fiber_id(nfiber);
    std::vector <double> fiber_center_x(nfiber);
    std::vector <double> fiber_center_y(nfiber);

    for (size_t j = 0; j < nfiber; ++j) {
        fiber_id[j] = hw_->fiber_id[j];
        double cx = hw_->fiber_pos_xy_mm[fiber_id[j]].first;
        double cy = hw_->fiber_pos_xy_mm[fiber_id[j]].second;
        fiber_center_x[j] = cx;
        fiber_center_y[j] = cy;
    }

    // shared_ptr reference counting is not threadsafe.  Here we extract
    // a copy of the "raw" pointers needed inside the parallel region.

    Targets * pobjs = objs.get();
    Tiles * ptiles = tiles.get();
    TargetTree * ptree = tree.get();
    Hardware * phw = hw_.get();

    #pragma omp parallel default(none) shared(logger, pobjs, ptiles, ptree, phw, ntile, tile_radius, nfiber, fiber_id, fiber_center_x, fiber_center_y, fiber_patrol_mm)
    {
        // We re-use these thread-local vectors to reduce the number
        // of times we are realloc'ing memory for every fiber.
        std::vector <int64_t> nearby;
        std::vector <KdTreePoint> nearby_tree_points;
        std::vector <int64_t> nearby_fiber;
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
                }
            }
            KDtree <KdTreePoint> nearby_tree(nearby_tree_points, 2);

            size_t fibers_with_targets = 0;

            double fiber_pos[2];

            for (size_t j = 0; j < nfiber; ++j) {
                thread_data[tid][fiber_id[j]].resize(0);
                fiber_pos[0] = fiber_center_x[j];
                fiber_pos[1] = fiber_center_y[j];
                auto fiber_xy = std::make_pair(fiber_center_x[j],
                                               fiber_center_y[j]);

                // Lookup targets near this fiber in focalplane
                // coordinates.
                nearby_fiber = nearby_tree.near(fiber_pos, 0.0,
                                                fiber_patrol_mm);

                if (nearby_fiber.size() == 0) {
                    // No targets for this fiber
                    continue;
                }

                // The kdtree gets us the targets that are close to our
                // region of interest.  Now go through these targets and
                // compute the precise distance from the fiber center.
                // We sort the available targets by priority now,
                // to avoid doing it multiple times later.
                result.clear();
                result_weight.clear();

                size_t tindx = 0;
                for (auto const & tnear : nearby_fiber) {
                    auto & obj = pobjs->data[tnear];
                    auto obj_xy = phw->radec2xy(tra, tdec, obj.ra,
                                                obj.dec);
                    double dist = geom::sq(fiber_xy, obj_xy);
                    if (dist > geom::sq(fiber_patrol_mm)) {
                        // outside the patrol radius
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
                            << ", fiber " << fiber_id[j]
                            << " append ID "
                            << result[wt.second] << " (type="
                            << (int)(pobjs->data[result[wt.second]].type) << ")"
                            << ", total priority " << wt.first;
                        logger.debug_tfg(tid, fiber_id[j],
                                         result[wt.second],
                                         logmsg.str().c_str());
                    }
                    thread_data[tid][fiber_id[j]].push_back(
                        result[wt.second]);
                }

                if (thread_data[tid][fiber_id[j]].size() > 0) {
                    fibers_with_targets++;
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

                for (size_t f = 0; f < nfiber; ++f) {
                    int32_t fid = fiber_id[f];
                    if (fmap.count(fid) == 0) {
                        // This fiber had no targets.
                        continue;
                    }
                    // logmsg.str("");
                    // logmsg << "tile " << ttile << ", fiber "
                    //     << fid << " copying thread result to output";
                    // logger.debug_tfg(ttile, fid, -1,
                    //                  logmsg.str().c_str());
                    data[ttile][fid] = fmap.at(fid);
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
        for (size_t j = 0; j < nfiber; ++j) {
            if (data.at(tid).count(fiber_id[j]) > 0) {
                total_avail += data.at(tid).at(fiber_id[j]).size();
            }
        }
        std::ostringstream msg;
        msg.str("");
        msg << "targets avail:  tile " << tid
            << ", " << total_avail << " total available targets";
        logger.debug(msg.str().c_str());
    }

    tm.stop();
    tm.report("Computing targets available to all tile / fibers");
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


fba::FibersAvailable::FibersAvailable(fba::TargetsAvailable::pshr tgsavail) {
    fba::Timer tm;
    tm.start();

    fba::Logger & logger = fba::Logger::get();

    data.clear();

    //std::cout << "FibersAvailable:  tgsavail has " << avail.size() << " tiles" << std::endl;

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
            for (auto const & favail : avail.at(tid)) {
                int32_t fbr = favail.first;
                for (auto const & tg : favail.second) {
                    if (thread_data.count(tg) == 0) {
                        thread_data[tg].resize(0);
                    }
                    if (logger.extra_debug()) {
                        logmsg.str("");
                        logmsg << "target " << tg
                            << " has available tile/fiber "
                            << tid << ", " << fbr;
                        logger.debug_tfg(tid, fbr, tg, logmsg.str().c_str());
                    }
                    thread_data[tg].push_back(std::make_pair(tid, fbr));
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
                    // logmsg.str("");
                    // logmsg << "target " << tg << " copy tile " << ft.first
                    //     << ", fiber " << ft.second
                    //     << " from thread-local memory";
                    // logger.debug_tfg(ft.first, ft.second, tg,
                    //                  logmsg.str().c_str());
                    data[tg].push_back(ft);
                }
            }
            thread_data.clear();
        }
    }

    tm.stop();
    tm.report("Computing tile / fibers available to all objects");
}


std::vector <std::pair <int32_t, int32_t> >
    fba::FibersAvailable::target_data(int64_t target) const {
    if (data.count(target) == 0) {
        return std::vector <std::pair <int32_t, int32_t> >();
    } else {
        return data.at(target);
    }
}
