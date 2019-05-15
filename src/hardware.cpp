// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <hardware.h>

#include <cmath>
#include <iostream>
#include <algorithm>
#include <sstream>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace fba = fiberassign;

namespace fbg = fiberassign::geom;


fba::Hardware::Hardware(std::vector <int32_t> const & location,
                        std::vector <int32_t> const & petal,
                        std::vector <int32_t> const & device,
                        std::vector <int32_t> const & slitblock,
                        std::vector <int32_t> const & blockfiber,
                        std::vector <int32_t> const & spectro,
                        std::vector <int32_t> const & fiber,
                        std::vector <int32_t> const & slit,
                        std::vector <std::string> const & device_type,
                        std::vector <double> const & x_mm,
                        std::vector <double> const & y_mm,
                        std::vector <double> const & z_mm,
                        std::vector <double> const & q_deg,
                        std::vector <double> const & s_mm,
                        std::vector <int32_t> const & status) {
    nloc = location.size();

    int32_t maxpetal = 0;
    for (auto const & p : petal) {
        if (p > maxpetal) {
            maxpetal = p;
        }
    }
    npetal = static_cast <size_t> (maxpetal + 1);

    loc_pos_xy_mm.clear();
    loc_pos_z_mm.clear();
    loc_pos_q_deg.clear();
    loc_pos_s_mm.clear();
    loc_petal.clear();
    loc_spectro.clear();
    loc_device.clear();
    loc_device_type.clear();
    loc_fiber.clear();
    loc_slit.clear();
    loc_slitblock.clear();
    loc_blockfiber.clear();
    locations.resize(0);
    state.clear();

    petal_locations.clear();
    for (int32_t i = 0; i < npetal; ++i) {
        petal_locations[i].resize(0);
    }

    for (int32_t i = 0; i < nloc; ++i) {
        locations.push_back(location[i]);
        loc_petal[location[i]] = petal[i];
        loc_spectro[location[i]] = spectro[i];
        loc_device[location[i]] = device[i];
        loc_device_type[location[i]] = device_type[i];
        loc_fiber[location[i]] = fiber[i];
        loc_slit[location[i]] = slit[i];
        loc_slitblock[location[i]] = slitblock[i];
        loc_blockfiber[location[i]] = blockfiber[i];
        petal_locations[petal[i]].push_back(location[i]);
        loc_pos_xy_mm[location[i]] = std::make_pair(x_mm[i], y_mm[i]);
        loc_pos_z_mm[location[i]] = z_mm[i];
        loc_pos_q_deg[location[i]] = q_deg[i];
        loc_pos_s_mm[location[i]] = s_mm[i];
        state[location[i]] = status[i];
        neighbors[location[i]].clear();
    }

    // Sort the locations
    std::stable_sort(locations.begin(), locations.end());
    for (int32_t i = 0; i < npetal; ++i) {
        std::stable_sort(petal_locations[i].begin(), petal_locations[i].end());
    }

    // hardcode for now...
    nfiber_petal = 500;

    // hardcode these for now...

    focalplane_radius_deg = 1.65;

    patrol_mm = 5.8;

    collide_mm = 1.98;

    collide_avg_mm = 3.2;

    no_collide_mm = 7.0;

    neighbor_radius_mm = 14.05;

    // FIXME:  this hardcoded number is taken from the original code.  Should
    // be verified...
    positioner_range = std::make_pair(3.0, 3.0);
    pos_theta_max_deg = 380.0;
    pos_phi_max_deg = 200.0;
    pos_theta_max = pos_theta_max_deg * M_PI / 180.0;
    pos_phi_max = pos_phi_max_deg * M_PI / 180.0;

    // Compute neighboring locations
    for (int32_t x = 0; x < nloc; ++x) {
        int32_t xid = locations[x];
        for (int32_t y = x + 1; y < nloc; ++y) {
            int32_t yid = locations[y];
            double dist = fbg::dist(loc_pos_xy_mm[xid],
                                    loc_pos_xy_mm[yid]);
            if (dist <= neighbor_radius_mm) {
                neighbors[xid].push_back(yid);
                neighbors[yid].push_back(xid);
            }
        }
    }

    ferrule_holder = create_ferrule_holder();
    central_body = create_central_body();
}


std::vector <int32_t> fba::Hardware::device_locations(
        std::string const & type) const {
    std::vector <int32_t> ret;
    for (auto const & fid : locations) {
        if (type.compare(loc_device_type.at(fid)) == 0) {
            ret.push_back(fid);
        }
    }
    return ret;
}


// Returns the radial distance on the focalplane (mm) given the angle,
// theta (radians).  This is simply a fit to the data provided.
double fba::Hardware::radial_ang2dist (double const & theta_rad) const {
    const double p[4] = {8.297e5, -1750., 1.394e4, 0.0};
    double dist_mm = 0.0;
    for (size_t i = 0; i < 4; ++i) {
        dist_mm = theta_rad * dist_mm + p[i];
    }
    return dist_mm;
}


// Returns the radial angle (theta) on the focalplane given the distance (mm)
double fba::Hardware::radial_dist2ang (double const & dist_mm) const {
    double delta_theta = 1e-4;
    double inv_delta = 1.0 / delta_theta;

    // starting guess
    double theta_rad = 0.01;

    double distcur;
    double distdelta;
    double correction;
    double error = 1.0;

    while (::abs(error) > 1e-7) {
        distcur = radial_ang2dist(theta_rad);
        distdelta = radial_ang2dist(theta_rad + delta_theta);
        error = distcur - dist_mm;
        correction = error / (inv_delta * (distdelta - distcur));
        theta_rad -= correction;
    }
    return theta_rad;
}


fiberassign::geom::dpair fba::Hardware::radec2xy(
    double const & tilera, double const & tiledec, double const & ra,
    double const & dec) const {

    double deg_to_rad = M_PI / 180.0;

    // Inclination is 90 degrees minus the declination in degrees
    double inc_rad = (90.0 - dec) * deg_to_rad;

    double ra_rad = ra * deg_to_rad;
    //double dec_rad = dec * deg_to_rad;
    double tilera_rad = tilera * deg_to_rad;
    double tiledec_rad = tiledec * deg_to_rad;

    double sin_inc_rad = ::sin(inc_rad);
    double x0 = sin_inc_rad * ::cos(ra_rad);
    double y0 = sin_inc_rad * ::sin(ra_rad);
    double z0 = ::cos(inc_rad);

    double cos_tilera_rad = ::cos(tilera_rad);
    double sin_tilera_rad = ::sin(tilera_rad);
    double x1 = cos_tilera_rad * x0 + sin_tilera_rad * y0;
    double y1 = -sin_tilera_rad * x0 + cos_tilera_rad * y0;
    double z1 = z0;

    double cos_tiledec_rad = ::cos(tiledec_rad);
    double sin_tiledec_rad = ::sin(tiledec_rad);
    double x = cos_tiledec_rad * x1 + sin_tiledec_rad * z1;
    double y = y1;
    double z = -sin_tiledec_rad * x1 + cos_tiledec_rad * z1;

    double ra_ang_rad = ::atan2(y, x);
    if (ra_ang_rad < 0) {
        ra_ang_rad = 2.0 * M_PI + ra_ang_rad;
    }

    double dec_ang_rad = (M_PI / 2)
        - ::acos(z / ::sqrt((x*x) + (y*y) + (z*z)) );

    double radius_rad = 2 *
        ::asin( ::sqrt( ::pow( ::sin(dec_ang_rad / 2), 2) +
        ::cos(dec_ang_rad) * ::pow( ::sin(ra_ang_rad / 2), 2) ) );

    double q_rad = ::atan2(z, -y);

    double radius_mm = radial_ang2dist(radius_rad);
    //std::cerr << "RADEC2XY:  radius_mm = " << radius_mm << std::endl;

    double x_focalplane = radius_mm * ::cos(q_rad);
    double y_focalplane = radius_mm * ::sin(q_rad);

    return std::make_pair(x_focalplane, y_focalplane);
}


void fba::Hardware::radec2xy_multi(
    double const & tilera, double const & tiledec,
    std::vector <double> const & ra,
    std::vector <double> const & dec,
    std::vector <std::pair <double, double> > & xy, int threads) const {

    size_t ntg = ra.size();
    xy.resize(ntg);

    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_num_threads();
    #endif
    int run_threads;
    if (threads > 0) {
        run_threads = threads;
    } else {
        run_threads = max_threads;
    }
    if (run_threads > max_threads) {
        run_threads = max_threads;
    }

    #pragma omp parallel for schedule(static) default(none) shared(ntg, tilera, tiledec, xy, ra, dec) num_threads(run_threads)
    for (size_t t = 0; t < ntg; ++t) {
        xy[t] = radec2xy(tilera, tiledec, ra[t], dec[t]);
    }

    return;
}


fiberassign::geom::dpair fba::Hardware::xy2radec(
        double const & tilera, double const & tiledec,
        double const & x_mm, double const & y_mm) const {

    double deg_to_rad = M_PI / 180.0;
    double rad_to_deg = 180.0 / M_PI;

    double tilera_rad = tilera * deg_to_rad;
    double tiledec_rad = tiledec * deg_to_rad;

    // radial distance on the focal plane
    double radius_mm = ::sqrt(x_mm * x_mm + y_mm * y_mm);
    double radius_rad = radial_dist2ang(radius_mm);

    // q is the angle the position makes with the +x-axis of focal plane
    double q_rad = ::atan2(y_mm, x_mm);
    //double q_deg = q_rad * 180.0 / M_PI;

    // The focal plane is oriented with +yfocal = +dec but +xfocal = -RA
    // Rotate clockwise around z by r_rad

    double x1 = ::cos(radius_rad);     // y0=0 so drop sin(r_rad) term
    double y1 = -::sin(radius_rad);    // y0=0 so drop cos(r_rad) term
    //double z1 = 0.0;

    // clockwise rotation around the x-axis

    double x2 = x1;
    double y2 = y1 * ::cos(q_rad);    // z1=0 so drop sin(q_rad) term
    double z2 = -y1 * ::sin(q_rad);

    double cos_tiledec = ::cos(tiledec_rad);
    double sin_tiledec = ::sin(tiledec_rad);
    double cos_tilera = ::cos(tilera_rad);
    double sin_tilera = ::sin(tilera_rad);

    // Clockwise rotation around y axis by declination of the tile center

    double x3 = cos_tiledec * x2 - sin_tiledec * z2;
    double y3 = y2;
    double z3 = sin_tiledec * x2 + cos_tiledec * z2;

    // Counter-clockwise rotation around the z-axis by the right
    // ascension of the tile center
    double x4 = cos_tilera * x3 - sin_tilera * y3;
    double y4 = sin_tilera * x3 + cos_tilera * y3;
    double z4 = z3;

    double ra_rad = ::atan2(y4, x4);
    if (ra_rad < 0.0) {
        ra_rad = 2.0 * M_PI + ra_rad;
    }

    double dec_rad = M_PI_2 - ::acos(z4);

    double ra = ::fmod( (ra_rad * rad_to_deg), 360.0);
    double dec = dec_rad * rad_to_deg;

    return std::make_pair(ra, dec);
}


void fba::Hardware::xy2radec_multi(
        double const & tilera, double const & tiledec,
        std::vector <double> const & x_mm, std::vector <double> const & y_mm,
        std::vector <std::pair <double, double> > & radec, int threads) const {
    size_t npos = x_mm.size();
    radec.resize(npos);

    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_num_threads();
    #endif
    int run_threads;
    if (threads > 0) {
        run_threads = threads;
    } else {
        run_threads = max_threads;
    }
    if (run_threads > max_threads) {
        run_threads = max_threads;
    }

    #pragma omp parallel for schedule(static) default(none) shared(npos, tilera, tiledec, x_mm, y_mm, radec) num_threads(run_threads)
    for (size_t i = 0; i < npos; ++i) {
        radec[i] = xy2radec(tilera, tiledec, x_mm[i], y_mm[i]);
    }

    return;
}


fbg::shape fba::Hardware::create_ferrule_holder() const {
    // Head disk
    fbg::circle head(std::make_pair(0.0, 0.0), 0.967);

    // Upper segment
    std::vector <fbg::dpair> points;
    points.push_back(std::make_pair(-3.0, 0.990));
    points.push_back(std::make_pair(0, 0.990));
    fbg::segments upper(points);

    // Lower segment
    points.clear();
    points.push_back(std::make_pair(-2.944, -1.339));
    points.push_back(std::make_pair(-2.944, -2.015));
    points.push_back(std::make_pair(-1.981, -1.757));
    points.push_back(std::make_pair(-1.844, -0.990));
    points.push_back(std::make_pair(0.0, -0.990));
    fbg::segments lower(points);

    // Construct the total shape
    fbg::circle_list circs;
    circs.push_back(head);
    fbg::segments_list segs;
    segs.push_back(upper);
    segs.push_back(lower);
    fbg::shape fh(std::make_pair(-3.0, 0.0), circs, segs);

    // Translate to proper location
    fh.transl(std::make_pair(6.0, 0.0));

    return fh;
}


fbg::shape fba::Hardware::create_central_body() const {
    // Disk
    fbg::circle disk(std::make_pair(3.0, 0.0), 2.095);

    // Segments
    std::vector <fbg::dpair> points;
    points.push_back(std::make_pair(5.095, -0.474));
    points.push_back(std::make_pair(4.358, -2.5));
    points.push_back(std::make_pair(2.771, -2.5));
    points.push_back(std::make_pair(1.759, -2.792));
    points.push_back(std::make_pair(0.905, -0.356));
    fbg::segments seg(points);

    // Total shape
    fbg::circle_list circs;
    circs.push_back(disk);
    fbg::segments_list segs;
    segs.push_back(seg);
    fbg::shape cb(std::make_pair(3.0, 0.0), circs, segs);

    return cb;
}


std::array <double, 4> fba::Hardware::pos_angles(
        fbg::dpair const & A, fbg::dpair const & posp) const {
    fba::Logger & logger = fba::Logger::get();
    std::array <double, 4> ang;
    double na = fbg::norm(A);
    double phi = 2.0 * ::acos(na / (2.0 * posp.first) );
    if ((phi < 0.0) || (phi > pos_phi_max)) {
        std::ostringstream o;
        o << "cannot move positioner to phi angle " << (phi * 180.0 / M_PI);
        logger.error(o.str().c_str());
    }
    // NOTE:  This was how theta was calculated in the original code.  The
    // reduction of theta range does not seem to make sense...
    // double theta = ::atan2(A.second, A.first) - (0.5 * phi)
    //     + (A.first < 0.0 ? M_PI : 0.0);
    double theta = ::atan2(A.second, A.first) - (0.5 * phi);
    // NOTE: theta range is currently 0-380 degrees, and atan2 return value
    // is in range +- 2PI.  So should not need a range check here...
    ang[0] = ::cos(theta);
    ang[1] = ::sin(theta);
    ang[2] = ::cos(phi);
    ang[3] = ::sin(phi);
    return ang;
}


void fba::Hardware::reposition_fiber(fbg::shape & cb, fbg::shape & fh,
                                     fbg::dpair const & center,
                                     fbg::dpair const & position,
                                     fbg::dpair const & posp) const {
    fba::Logger & logger = fba::Logger::get();
    // repositions positioner (central body, ferule holder)
    fbg::dpair offset = std::make_pair(position.first - center.first,
                                       position.second - center.second);
    if (fbg::sq(posp.first + posp.second) < fbg::sq(offset)) {
        std::ostringstream o;
        o << "cannot move fiber at (" << center.first << ", "
            << center.second << ") to position (" << position.first
            << ", " << position.second << ")";
        logger.error(o.str().c_str());
    }

    std::array <double, 4> ang = pos_angles(offset, posp);
    fbg::dpair theta = std::make_pair(ang[0], ang[1]);
    if (std::isnan(theta.first) || std::isnan(theta.second)) {
        std::ostringstream o;
        o << "Moving fiber at (" << center.first << ", "
            << center.second << ") to position (" << position.first
            << ", " << position.second << "): theta angle is NaN";
        logger.error(o.str().c_str());
    }
    fbg::dpair phi = std::make_pair(ang[2], ang[3]);
    if (std::isnan(phi.first) || std::isnan(phi.second)) {
        std::ostringstream o;
        o << "Moving fiber at (" << center.first << ", "
            << center.second << ") to position (" << position.first
            << ", " << position.second << "): phi angle is NaN";
        logger.error(o.str().c_str());
    }

    cb.rotation_origin(theta);
    fh.rotation_origin(theta);
    fh.rotation(phi);
    fh.transl(center);
    cb.transl(center);
    return;
}


bool fba::Hardware::collide(fbg::dpair center1, fbg::dpair position1,
                            fbg::dpair center2, fbg::dpair position2) const {
    // Check for collisions if we move fibers at given central positions
    // (in mm) to the x/y positions specified.
    double dist = fbg::dist(position1, position2);

    if (dist < collide_mm) {
        // Guaranteed to collide.
        return true;
    }

    if (dist > no_collide_mm) {
        // Guaranteed to NOT collide.
        return false;
    }

    // We need to do a detailed calculation of the fiber holder and central
    // body positions and test for intersection.

    fbg::shape fh1(ferrule_holder);
    fbg::shape fh2(ferrule_holder);
    fbg::shape cb1(central_body);
    fbg::shape cb2(central_body);

    reposition_fiber(cb1, fh1, center1, position1, positioner_range);
    reposition_fiber(cb2, fh2, center2, position2, positioner_range);

    if (fbg::intersect(fh1, fh2)) {
        return true;
    }
    if (fbg::intersect(cb1, fh2)) {
        return true;
    }
    if (fbg::intersect(cb2, fh1)) {
        return true;
    }

    return false;
}


std::pair <fbg::shape, fbg::shape> fba::Hardware::loc_position(
    int32_t loc, fbg::dpair const & xy) const {

    // Construct a positioner shape for a given fiber and location.
    fbg::shape fh(ferrule_holder);
    fbg::shape cb(central_body);

    auto const & center = loc_pos_xy_mm.at(loc);

    reposition_fiber(cb, fh, center, xy, positioner_range);
    return std::make_pair(cb, fh);
}


std::vector <std::pair <fbg::shape, fbg::shape> >
    fba::Hardware::loc_position_multi(
        std::vector <int32_t> const & loc,
        std::vector <fbg::dpair> const & xy, int threads) const {
    size_t nlc = loc.size();
    std::vector <std::pair <fbg::shape, fbg::shape> > result(nlc);

    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_num_threads();
    #endif
    int run_threads = 1;
    if (threads > 0) {
        run_threads = threads;
    } else {
        run_threads = max_threads;
    }
    if (run_threads > max_threads) {
        run_threads = max_threads;
    }

    #pragma omp parallel for schedule(static) default(none) shared(nlc, loc, xy, result) num_threads(run_threads)
    for (size_t f = 0; f < nlc; ++f) {
        int32_t fid = loc[f];
        // Construct a positioner shape for a given fiber and location.
        fbg::shape fh(ferrule_holder);
        fbg::shape cb(central_body);
        auto const & center = loc_pos_xy_mm.at(fid);
        reposition_fiber(cb, fh, center, xy[f], positioner_range);
        result[f] = std::make_pair(cb, fh);
    }
    return result;
}


std::vector <bool> fba::Hardware::check_collisions_xy(
    std::vector <int32_t> const & loc,
    std::vector <fbg::dpair> const & xy, int threads) const {

    size_t nlc = loc.size();

    auto fpos = loc_position_multi(loc, xy);

    std::map <int32_t, int32_t> loc_indx;

    // Build list of all location pairs to check for a collision
    std::map <int32_t, std::vector <int32_t> > checklookup;
    size_t idx = 0;
    for (auto const & lid : loc) {
        loc_indx[lid] = idx;
        for (auto const & nb : neighbors.at(lid)) {
            int32_t low;
            int32_t high;
            if (lid < nb) {
                low = lid;
                high = nb;
            } else {
                low = nb;
                high = lid;
            }
            if (checklookup.count(low) == 0) {
                checklookup[low].clear();
            }
            bool found = false;
            for (auto const & ck : checklookup[low]) {
                if (ck == high) {
                    found = true;
                }
            }
            if (! found) {
                checklookup[low].push_back(high);
            }
        }
        idx++;
    }
    std::vector <std::pair <int32_t, int32_t> > checkpairs;
    for (auto const & it : checklookup) {
        int32_t low = it.first;
        for (auto const & high : it.second) {
            checkpairs.push_back(std::make_pair(low, high));
        }
    }
    checklookup.clear();

    size_t npairs = checkpairs.size();

    std::vector <bool> result(nlc);
    result.assign(nlc, false);

    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_num_threads();
    #endif
    int run_threads = 1;
    if (threads > 0) {
        run_threads = threads;
    } else {
        run_threads = max_threads;
    }
    if (run_threads > max_threads) {
        run_threads = max_threads;
    }

    #pragma omp parallel for schedule(static) default(none) shared(npairs, checkpairs, fpos, result, loc_indx) num_threads(run_threads)
    for (size_t p = 0; p < npairs; ++p) {
        int32_t flow = checkpairs[p].first;
        int32_t fhigh = checkpairs[p].second;
        bool hit = false;
        auto const & cb1 = fpos[loc_indx[flow]].first;
        auto const & fh1 = fpos[loc_indx[flow]].second;
        auto const & cb2 = fpos[loc_indx[fhigh]].first;
        auto const & fh2 = fpos[loc_indx[fhigh]].second;
        if (fbg::intersect(fh1, fh2)) {
            hit = true;
        }
        if (fbg::intersect(cb1, fh2)) {
            hit = true;
        }
        if (fbg::intersect(cb2, fh1)) {
            hit = true;
        }
        if (hit) {
            #pragma omp critical
            {
                result[flow] = true;
                result[fhigh] = true;
            }
        }
    }
    return result;
}


std::vector <bool> fba::Hardware::check_collisions_thetaphi(
    std::vector <int32_t> const & loc,
    std::vector <double> const & theta,
    std::vector <double> const & phi, int threads) const {

    size_t nlc = loc.size();

    std::map <int32_t, int32_t> loc_indx;

    // Build list of all location pairs to check for a collision
    std::map <int32_t, std::vector <int32_t> > checklookup;
    size_t idx = 0;
    for (auto const & lid : loc) {
        loc_indx[lid] = idx;
        for (auto const & nb : neighbors.at(lid)) {
            int32_t low;
            int32_t high;
            if (lid < nb) {
                low = lid;
                high = nb;
            } else {
                low = nb;
                high = lid;
            }
            if (checklookup.count(low) == 0) {
                checklookup[low].clear();
            }
            bool found = false;
            for (auto const & ck : checklookup[low]) {
                if (ck == high) {
                    found = true;
                }
            }
            if (! found) {
                checklookup[low].push_back(high);
            }
        }
        idx++;
    }
    std::vector <std::pair <int32_t, int32_t> > checkpairs;
    for (auto const & it : checklookup) {
        int32_t low = it.first;
        for (auto const & high : it.second) {
            checkpairs.push_back(std::make_pair(low, high));
        }
    }
    checklookup.clear();

    size_t npairs = checkpairs.size();

    std::vector <bool> result(nlc);
    result.assign(nlc, false);

    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_num_threads();
    #endif
    int run_threads = 1;
    if (threads > 0) {
        run_threads = threads;
    } else {
        run_threads = max_threads;
    }
    if (run_threads > max_threads) {
        run_threads = max_threads;
    }

    #pragma omp parallel for schedule(static) default(none) shared(npairs, checkpairs, result, loc_indx, theta, phi) num_threads(run_threads)
    for (size_t p = 0; p < npairs; ++p) {
        int32_t flow = checkpairs[p].first;
        int32_t fhigh = checkpairs[p].second;
        double thetalow = theta[loc_indx[flow]];
        double philow = phi[loc_indx[flow]];
        double thetahigh = theta[loc_indx[fhigh]];
        double phihigh = phi[loc_indx[fhigh]];

        fbg::shape fhlow(ferrule_holder);
        fbg::shape fhhigh(ferrule_holder);
        fbg::shape cblow(central_body);
        fbg::shape cbhigh(central_body);

        double costheta = ::cos(thetalow);
        double sintheta = ::sin(thetalow);
        double cosphi = ::cos(philow);
        double sinphi = ::sin(philow);
        std::pair <double, double> cossintheta = std::make_pair(costheta,
                                                                sintheta);
        std::pair <double, double> cossinphi = std::make_pair(cosphi,
                                                              sinphi);
        cblow.rotation_origin(cossintheta);
        fhlow.rotation_origin(cossintheta);
        fhlow.rotation(cossinphi);
        fhlow.transl(loc_pos_xy_mm.at(flow));
        cblow.transl(loc_pos_xy_mm.at(flow));

        costheta = ::cos(thetahigh);
        sintheta = ::sin(thetahigh);
        cosphi = ::cos(phihigh);
        sinphi = ::sin(phihigh);
        cossintheta = std::make_pair(costheta, sintheta);
        cossinphi = std::make_pair(cosphi, sinphi);
        cbhigh.rotation_origin(cossintheta);
        fhhigh.rotation_origin(cossintheta);
        fhhigh.rotation(cossinphi);
        fhhigh.transl(loc_pos_xy_mm.at(fhigh));
        cbhigh.transl(loc_pos_xy_mm.at(fhigh));

        bool hit = false;
        if (fbg::intersect(fhlow, fhhigh)) {
            hit = true;
        }
        if (fbg::intersect(cblow, fhhigh)) {
            hit = true;
        }
        if (fbg::intersect(cblow, fhhigh)) {
            hit = true;
        }
        if (hit) {
            #pragma omp critical
            {
                result[flow] = true;
                result[fhigh] = true;
            }
        }
    }
    return result;
}
