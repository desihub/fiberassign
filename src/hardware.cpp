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


fba::Hardware::Hardware(std::string const & timestr,
                        std::vector <int32_t> const & location,
                        std::vector <int32_t> const & petal,
                        std::vector <int32_t> const & device,
                        std::vector <int32_t> const & slitblock,
                        std::vector <int32_t> const & blockfiber,
                        std::vector <int32_t> const & fiber,
                        std::vector <std::string> const & device_type,
                        std::vector <double> const & x_mm,
                        std::vector <double> const & y_mm,
                        std::vector <int32_t> const & status,
                        std::vector <double> const & theta_offset,
                        std::vector <double> const & theta_min,
                        std::vector <double> const & theta_max,
                        std::vector <double> const & theta_pos,
                        std::vector <double> const & theta_arm,
                        std::vector <double> const & phi_offset,
                        std::vector <double> const & phi_min,
                        std::vector <double> const & phi_max,
                        std::vector <double> const & phi_pos,
                        std::vector <double> const & phi_arm,
                        std::vector <double> const & ps_radius,
                        std::vector <double> const & ps_theta,
                        std::vector <double> const & arclen,
                        std::vector <fbg::shape> const & excl_theta,
                        std::vector <fbg::shape> const & excl_phi,
                        std::vector <fbg::shape> const & excl_gfa,
                        std::vector <fbg::shape> const & excl_petal) {
    nloc = location.size();

    fba::Logger & logger = fba::Logger::get();
    std::ostringstream logmsg;

    timestr_ = timestr;

    ps_radius_ = ps_radius;
    ps_theta_ = ps_theta;
    ps_size_ = ps_radius.size();
    arclen_ = arclen;

    int32_t maxpetal = 0;
    for (auto const & p : petal) {
        if (p > maxpetal) {
            maxpetal = p;
        }
    }
    npetal = static_cast <size_t> (maxpetal + 1);

    int32_t maxslitblock = 0;
    for (auto const & p : slitblock) {
        if (p > maxslitblock) {
            maxslitblock = p;
        }
    }
    nslitblock = static_cast <size_t> (maxslitblock + 1);

    loc_pos_cs5_mm.clear();
    loc_pos_curved_mm.clear();
    loc_petal.clear();
    loc_device.clear();
    loc_device_type.clear();
    loc_fiber.clear();
    loc_slitblock.clear();
    loc_blockfiber.clear();
    locations.resize(0);
    state.clear();
    loc_theta_offset.clear();
    loc_theta_min.clear();
    loc_theta_max.clear();
    loc_theta_pos.clear();
    loc_theta_arm.clear();
    loc_phi_offset.clear();
    loc_phi_min.clear();
    loc_phi_max.clear();
    loc_phi_pos.clear();
    loc_phi_arm.clear();
    loc_theta_excl.clear();
    loc_phi_excl.clear();
    petal_edge.clear();
    gfa_edge.clear();

    petal_locations.clear();
    for (int32_t i = 0; i < npetal; ++i) {
        petal_locations[i].resize(0);
    }

    int32_t stcount = 0;

    for (int32_t i = 0; i < nloc; ++i) {
        int32_t loc = location[i];
        locations.push_back(loc);
        loc_petal[loc] = petal[i];
        loc_device[loc] = device[i];
        loc_device_type[loc] = device_type[i];
        loc_fiber[loc] = fiber[i];
        loc_slitblock[loc] = slitblock[i];
        loc_blockfiber[loc] = blockfiber[i];
        petal_locations[petal[i]].push_back(loc);
        loc_pos_cs5_mm[loc] = std::make_pair(x_mm[i], y_mm[i]);
        state[loc] = status[i];
        if ((status[i] & FIBER_STATE_STUCK) || (status[i] & FIBER_STATE_BROKEN)) {
            stcount++;
        }
        neighbors[loc].clear();
        petal_edge[loc] = false;
        gfa_edge[loc] = false;
        loc_theta_offset[loc] = theta_offset[i] * M_PI / 180.0;
        loc_theta_min[loc] = theta_min[i] * M_PI / 180.0;
        loc_theta_max[loc] = theta_max[i] * M_PI / 180.0;
        loc_theta_pos[loc] = theta_pos[i] * M_PI / 180.0;
        loc_theta_arm[loc] = theta_arm[i];
        loc_phi_offset[loc] = phi_offset[i] * M_PI / 180.0;
        loc_phi_min[loc] = phi_min[i] * M_PI / 180.0;
        loc_phi_max[loc] = phi_max[i] * M_PI / 180.0;
        loc_phi_pos[loc] = phi_pos[i] * M_PI / 180.0;
        loc_phi_arm[loc] = phi_arm[i];

        // For stuck positioners or broken fibers, the current fixed theta / phi angle
        // values should be specified.  Other movable positioners (even if restricted)
        // will have these values set to zero and they will not be used in fiberassign.
        // Check that these fixed angles are within range and print a message if not.
        // This can happen when a positioner is within its "physical" angle range but
        // not within the "targetable" range we get from desimodel.  Internally for
        // these fixed positioners we just take the given positions as true, regardless
        // of range.
        if (state[loc] & (FIBER_STATE_STUCK | FIBER_STATE_BROKEN)) {
            if ((loc_theta_pos[loc] < loc_theta_min[loc]) ||
                (loc_theta_pos[loc] > loc_theta_max[loc])) {
                logmsg.str("");
                logmsg << ((state[loc] & FIBER_STATE_STUCK) ? "STUCK" : "     ")
                       << " | "
                       << ((state[loc] & FIBER_STATE_BROKEN) ? "BROKEN" : "      ")
                       << " positioner (loc " << loc << ") theta value is outside of range: theta = "
                       << loc_theta_pos[loc] << " vs range [" << loc_theta_min[loc]
                       << ", " << loc_theta_max[loc] << "].";
                logger.debug(logmsg.str().c_str());
            }
            if ((loc_phi_pos[loc] < loc_phi_min[loc]) ||
                (loc_phi_pos[loc] > loc_phi_max[loc])) {
                logmsg.str("");
                logmsg << ((state[loc] & FIBER_STATE_STUCK) ? "STUCK" : "     ")
                       << " | "
                       << ((state[loc] & FIBER_STATE_BROKEN) ? "BROKEN" : "      ")
                       << " positioner (loc " << loc << ") phi value is outside of range: phi = "
                       << loc_phi_pos[loc] << " vs range [" << loc_phi_min[loc]
                       << ", " << loc_phi_max[loc] << "].";
                logger.debug(logmsg.str().c_str());
            }
        }

        loc_theta_excl[loc] = excl_theta[i];
        loc_phi_excl[loc] = excl_phi[i];
        loc_gfa_excl[loc] = excl_gfa[i];
        loc_petal_excl[loc] = excl_petal[i];
    }

    logmsg.str("");
    logmsg << "Focalplane has " << stcount << " fibers that are stuck or broken";
    logger.info(logmsg.str().c_str());

    // Sort the locations
    std::stable_sort(locations.begin(), locations.end());
    for (int32_t i = 0; i < npetal; ++i) {
        std::stable_sort(petal_locations[i].begin(), petal_locations[i].end());
    }

    // Hard-coded parameters.  These could be moved to desimodel and passed
    // into this constructor as arguments.

    // The number of science positioners per petal.
    nfiber_petal = 500;

    // The tile / focalplane radius in degrees, used for selecting targets
    // that are available to a particular tile.
    focalplane_radius_deg = 1.65;

    // The radius in mm on the focalplane for considering which positioners
    // are "neighbors".
    neighbor_radius_mm = 14.05;

    // The amount to reduce the total arm length when considering which
    // targets are reachable by a positioner.  This was set to 200 microns
    // long ago...
    patrol_buffer_mm = 0.2;

    // Compute the positioner locations in the curved focal surface.
    for (int32_t i = 0; i < nloc; ++i) {
        int32_t lid = locations[i];
        double xcs5 = loc_pos_cs5_mm.at(lid).first;
        double ycs5 = loc_pos_cs5_mm.at(lid).second;
        double rot = ::atan2(ycs5, xcs5);
        double cs5_mm = ::sqrt(xcs5 * xcs5 + ycs5 * ycs5);
        double theta_rad = radial_dist2ang_CS5(cs5_mm);
        double curve_mm = radial_ang2dist_curved(theta_rad);
        double xcurve = curve_mm * ::cos(rot);
        double ycurve = curve_mm * ::sin(rot);
        loc_pos_curved_mm[lid] = std::make_pair(xcurve, ycurve);
    }

    // Compute neighboring locations
    for (int32_t x = 0; x < nloc; ++x) {
        int32_t xid = locations[x];
        for (int32_t y = x + 1; y < nloc; ++y) {
            int32_t yid = locations[y];
            double dist = fbg::dist(loc_pos_curved_mm[xid],
                                    loc_pos_curved_mm[yid]);
            if (dist <= neighbor_radius_mm) {
                neighbors[xid].push_back(yid);
                neighbors[yid].push_back(xid);
            }
        }
    }

    // For each location, we rotate the petal and GFA exclusion polygons
    // to the correct petal location.

    for (auto const & lid : locations) {
        int32_t pt = loc_petal[lid];
        double petalrot_deg = fmod((double)(7 + pt) * 36.0, 360.0);
        double petalrot_rad = petalrot_deg * M_PI / 180.0;
        auto csang = std::make_pair(cos(petalrot_rad), sin(petalrot_rad));
        loc_gfa_excl.at(lid).rotation_origin(csang);
        loc_petal_excl.at(lid).rotation_origin(csang);
    }
}


std::string fba::Hardware::time() const {
    return timestr_;
}


std::vector <double> fba::Hardware::platescale_radius_mm() const {
    return ps_radius_;
}


std::vector <double> fba::Hardware::platescale_theta_deg() const {
    return ps_theta_;
}


std::vector <double> fba::Hardware::radial_arclen() const {
    return arclen_;
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

// Small helper function to seek to the correct elements for linear interpolation.
void helper_vec_seek(
    std::vector <double> const & data,
    double const & val,
    size_t & ilow,
    size_t & ihigh
) {
    std::vector <double>::const_iterator bound =
        std::lower_bound(data.begin(), data.end(), val);
    if (bound - data.end() == 0) {
        // We are extrapolating off the end
        ilow = data.size() - 2;
        ihigh = data.size() - 1;
    } else if (bound - data.begin() == 0) {
        // Extrapolating at the beginning
        ilow = 0;
        ihigh = 1;
    } else {
        // Somewhere in the middle
        ilow = bound - data.begin();
        if (data[ilow] > val) {
            ilow--;
        }
        ihigh = ilow + 1;
    }
    return;
}


// Returns the radial distance in CS5 on the focalplane (mm) given the angle,
// theta (radians).  This does a linear interpolation to the platescale
// data.
double fba::Hardware::radial_ang2dist_CS5 (double const & theta_rad) const {
    // platescale data is in degrees
    double theta_deg = theta_rad * 180.0 / M_PI;

    // Seek to correct entry
    size_t ilow;
    size_t ihigh;
    helper_vec_seek(ps_theta_, theta_deg, ilow, ihigh);

    // Interpolate
    double xfrac = (theta_deg - ps_theta_[ilow]) / (ps_theta_[ihigh] - ps_theta_[ilow]);
    double dist_mm = ps_radius_[ilow] + xfrac * (ps_radius_[ihigh] - ps_radius_[ilow]);

    return dist_mm;
}


// Returns the radial angle (theta) on the focalplane given the distance (mm) in CS5.
// This does a linear interpolation to the platescale data.
double fba::Hardware::radial_dist2ang_CS5 (double const & dist_mm) const {
    // Seek to correct entry
    size_t ilow;
    size_t ihigh;
    helper_vec_seek(ps_radius_, dist_mm, ilow, ihigh);

    // Interpolate
    double xfrac = (dist_mm - ps_radius_[ilow]) /
        (ps_radius_[ihigh] - ps_radius_[ilow]);
    double theta_deg = ps_theta_[ilow] + xfrac * (ps_theta_[ihigh] - ps_theta_[ilow]);

    // platescale data is in degrees
    double theta_rad = theta_deg * M_PI / 180.0;

    return theta_rad;
}


// Returns the radial arc length S(R) on the focal surface (mm) given the angle,
// theta (radians).  This does a linear interpolation to the model.
double fba::Hardware::radial_ang2dist_curved (double const & theta_rad) const {
    // platescale data is in degrees
    double theta_deg = theta_rad * 180.0 / M_PI;

    // Seek to correct entry
    size_t ilow;
    size_t ihigh;
    helper_vec_seek(ps_theta_, theta_deg, ilow, ihigh);

    // Interpolate
    double xfrac = (theta_deg - ps_theta_[ilow]) / (ps_theta_[ihigh] - ps_theta_[ilow]);
    double arc_mm = arclen_[ilow] + xfrac * (arclen_[ihigh] - arclen_[ilow]);

    return arc_mm;
}


// Returns the radial angle (theta) on the focal surface given the arc length (mm).
// This does a linear interpolation to the model.
double fba::Hardware::radial_dist2ang_curved (double const & arc_mm) const {
    // Seek to correct entry
    size_t ilow;
    size_t ihigh;
    helper_vec_seek(ps_radius_, arc_mm, ilow, ihigh);

    // Interpolate
    double xfrac = (arc_mm - arclen_[ilow]) /
        (arclen_[ihigh] - arclen_[ilow]);
    double theta_deg = ps_theta_[ilow] + xfrac * (ps_theta_[ihigh] - ps_theta_[ilow]);

    // platescale data is in degrees
    double theta_rad = theta_deg * M_PI / 180.0;

    return theta_rad;
}


fiberassign::geom::dpair fba::Hardware::radec2xy(
    double const & tilera, double const & tiledec, double const & tiletheta,
    double const & ra, double const & dec, bool use_CS5) const {

    double deg_to_rad = M_PI / 180.0;

    // Inclination is 90 degrees minus the declination in degrees
    double inc_rad = (90.0 - dec) * deg_to_rad;

    double ra_rad = ra * deg_to_rad;
    //double dec_rad = dec * deg_to_rad;
    double tilera_rad = tilera * deg_to_rad;
    double tiledec_rad = tiledec * deg_to_rad;
    double tiletheta_rad = tiletheta * deg_to_rad;

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

    double radius_mm;
    if (use_CS5) {
        radius_mm = radial_ang2dist_CS5(radius_rad);
    } else {
        radius_mm = radial_ang2dist_curved(radius_rad);
    }

    // Apply field rotation
    double rotated = q_rad + tiletheta_rad;

    double x_focalplane = radius_mm * ::cos(rotated);
    double y_focalplane = radius_mm * ::sin(rotated);

    return std::make_pair(x_focalplane, y_focalplane);
}


void fba::Hardware::radec2xy_multi(
    double const & tilera, double const & tiledec, double const & tiletheta,
    std::vector <double> const & ra,
    std::vector <double> const & dec,
    std::vector <std::pair <double, double> > & xy, bool use_CS5,
    int threads) const {

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

    #pragma omp parallel for schedule(static) default(none) shared(ntg, tilera, tiledec, tiletheta, xy, ra, dec, use_CS5) num_threads(run_threads)
    for (size_t t = 0; t < ntg; ++t) {
        xy[t] = radec2xy(tilera, tiledec, tiletheta, ra[t], dec[t], use_CS5);
    }

    return;
}


fiberassign::geom::dpair fba::Hardware::xy2radec(
        double const & tilera, double const & tiledec,
        double const & tiletheta,
        double const & x_mm, double const & y_mm, bool use_CS5) const {

    double deg_to_rad = M_PI / 180.0;
    double rad_to_deg = 180.0 / M_PI;

    double tilera_rad = tilera * deg_to_rad;
    double tiledec_rad = tiledec * deg_to_rad;
    double tiletheta_rad = tiletheta * deg_to_rad;

    // radial distance on the focal plane
    double radius_mm = ::sqrt(x_mm * x_mm + y_mm * y_mm);

    double radius_rad;
    if (use_CS5) {
        radius_rad = radial_dist2ang_CS5(radius_mm);
    } else {
        radius_rad = radial_dist2ang_curved(radius_mm);
    }

    // q is the angle the position makes with the +x-axis of focal plane
    double rotated = ::atan2(y_mm, x_mm);

    // Remove field rotation
    double q_rad = rotated - tiletheta_rad;

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
        double const & tiletheta,
        std::vector <double> const & x_mm, std::vector <double> const & y_mm,
        std::vector <std::pair <double, double> > & radec,
        bool use_CS5, int threads) const {
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

    #pragma omp parallel for schedule(static) default(none) shared(npos, tilera, tiledec, tiletheta, x_mm, y_mm, radec, use_CS5) num_threads(run_threads)
    for (size_t i = 0; i < npos; ++i) {
        radec[i] = xy2radec(tilera, tiledec, tiletheta, x_mm[i], y_mm[i], use_CS5);
    }

    return;
}


double _angle_diff(double hi, double low) {
    double twopi = 2.0 * M_PI;
    // range reduction to [-Pi, Pi)
    if (hi >= M_PI) {
        hi -= twopi;
    }
    if (hi < -M_PI) {
        hi += twopi;
    }
    if (low >= M_PI) {
        low -= twopi;
    }
    if (low < -M_PI) {
        low += twopi;
    }
    double diff = hi - low;
    if (diff > 0.0) {
        return diff;
    } else {
        return (twopi + diff);
    }
}


bool _outside_theta_phi_range(
    double const & theta, double const & phi,
    double const & theta_zero, double const & theta_min, double const & theta_max,
    double const & phi_zero, double const & phi_min, double const & phi_max
) {
    // The definitions of the theta / phi angles and ranges can be found in
    // DESI-0899.

    // Theta angle.  The theta offset loaded from desimodel is already in global
    // focalplane coordinates (not petal-local).  The theta min / max values are
    // relative to this offset (not the coordinate system).  Compute the angle
    // difference rotating both directions from theta_zero and see if we can reach
    // it at least one of those ways.  Note that theta_min is negative to indicate
    // a clockwise rotation from theta_zero.
    double diff_hi = _angle_diff(theta, theta_zero);
    double diff_lo = _angle_diff(theta_zero, theta);

    // Since calling code may add theta_offset (aka theta_zero here),
    // which we then subtract off, we want to put a little margin on
    // the min/max checks.  (In most cases where we are passing
    // theta/phi values through from the desimodel inputs and we need
    // to add offsets, we also ignore this range check, so this is
    // maybe overly cautious.)
    double eps = 1e-12;

    if ((   ( diff_hi > (theta_max + eps)) || ( diff_hi < (theta_min - eps)))
        && ((-diff_lo > (theta_max + eps)) || (-diff_lo < (theta_min - eps)))
        ) {
        return true;
    }

    // Phi angle.  The phi offset (phi_zero) is relative to the coordinate system
    // defined by the theta arm along the X axis.  The phi min / max values are
    // relative to this offset (not the theta-arm coordinate system).  Note that
    // a negative phi_min indicates a clockwise rotation from phi_zero.
    diff_hi = _angle_diff(phi, phi_zero);
    diff_lo = _angle_diff(phi_zero, phi);
    if ((   ( diff_hi > (phi_max + eps)) || ( diff_hi < (phi_min - eps)))
        && ((-diff_lo > (phi_max + eps)) || (-diff_lo < (phi_min - eps)))
        ) {
        return true;
    }

    // not outside range
    return false;
}


bool fba::Hardware::move_positioner_thetaphi(
        fbg::shape & shptheta, fbg::shape & shpphi,
        fbg::dpair const & center, double theta, double phi,
        double theta_arm, double phi_arm, double theta_zero, double phi_zero,
        double theta_min, double phi_min, double theta_max, double phi_max,
        bool ignore_range
        ) const {
    // Check that requested angles are in range.
    if (
        _outside_theta_phi_range(
            theta, phi, theta_zero, theta_min, theta_max, phi_zero, phi_min, phi_max
        )
    ) {
        if (!ignore_range)
            return true;
    }

    double ctheta = ::cos(theta);
    double stheta = ::sin(theta);
    double cphi = ::cos(phi);
    double sphi = ::sin(phi);
    auto cstheta = std::make_pair(ctheta, stheta);
    auto csphi = std::make_pair(cphi, sphi);

    // Move the phi polygon into the fully extended position along the X axis.
    shpphi.transl(std::make_pair(theta_arm, 0.0));

    // std::cout << "move_positioner_thetaphi:  after transl:" << std::endl;
    // shptheta.print();
    // shpphi.print();

    // Rotate fully extended positioner an angle of theta about the center.
    shptheta.rotation_origin(cstheta);
    shpphi.rotation_origin(cstheta);

    // std::cout << "move_positioner_thetaphi:  after rot origin:" << std::endl;
    // shptheta.print();
    // shpphi.print();

    // Rotate just the phi arm an angle phi about the theta arm center.
    shpphi.rotation(csphi);

    // std::cout << "move_positioner_thetaphi:  after phi rot of " << csphi.first << ", " << csphi.second << ":" << std::endl;
    // shptheta.print();
    // shpphi.print();

    // Translate the whole positioner to the center.
    shpphi.transl(center);
    shptheta.transl(center);

    // std::cout << "move_positioner_thetaphi:  after center transl:" << std::endl;
    // shptheta.print();
    // shpphi.print();

    return false;
}


bool fba::Hardware::xy_to_thetaphi(
        double & theta, double & phi,
        fbg::dpair const & center, fbg::dpair const & position,
        double theta_arm, double phi_arm, double theta_zero, double phi_zero,
        double theta_min, double phi_min, double theta_max, double phi_max
        ) const {
    fbg::dpair offset = std::make_pair(position.first - center.first,
                                       position.second - center.second);

    double sq_theta_arm = theta_arm * theta_arm;
    double sq_phi_arm = phi_arm * phi_arm;
    double sq_offset = fbg::sq(offset);
    double sq_total_arm = fbg::sq(theta_arm + phi_arm);
    double sq_diff_arm = fbg::sq(theta_arm - phi_arm);

    phi = M_PI;
    theta = 0.0;

    if (fabs(sq_offset - sq_total_arm) <=
        std::numeric_limits<float>::epsilon()) {
        // We are at the maximum arm extension.  Force phi angle to zero
        // and compute theta.
        phi = 0.0;
        theta = ::atan2(offset.second, offset.first);
    } else if (fabs(sq_diff_arm - sq_offset) <=
        std::numeric_limits<float>::epsilon()) {
        // We are at the limit of the arm folded inwards.  Force phi angle
        // to PI and compute theta.
        phi = M_PI;
        theta = ::atan2(offset.second, offset.first);
    } else {
        // We are on neither limit.

        if (sq_total_arm < sq_offset) {
            // Physically impossible to get there for any choice of angles
            // std::cout << "xy_to_thetaphi: sqoffset - sqtot = " << sq_offset - sq_total_arm << std::endl;
            return true;
        }

        if (sq_offset < sq_diff_arm) {
            // Physically impossible to get there for any choice of angles
            // std::cout << "xy_to_thetaphi: sqdiff - sqoffset = " << sq_diff_arm - sq_offset << std::endl;
            return true;
        }

        // Use law of cosines to compute "opening" angle at the "elbow".
        double opening = ::acos((sq_theta_arm + sq_phi_arm - sq_offset)
                                / (2.0 * theta_arm * phi_arm));

        // std::cout << "xy_to_thetaphi: opening = " << opening << std::endl;

        // The PHI angle is just the supplement of this.
        phi = M_PI - opening;

        // std::cout << "xy_to_thetaphi: phi = " << phi << std::endl;

        // Compute the theta angle.
        // Use law of cosines to compute angle from theta arm to the line from
        // the origin to the X/Y position.
        double nrm_offset = ::sqrt(sq_offset);
        double txy = ::acos((sq_theta_arm + sq_offset - sq_phi_arm)
                            / (2 * theta_arm * nrm_offset));

        // std::cout << "xy_to_thetaphi: txy = " << txy << std::endl;

        theta = ::atan2(offset.second, offset.first) - txy;
    }
    // Check that requested angles are in range.
    if (
        _outside_theta_phi_range(
            theta, phi, theta_zero, theta_min, theta_max, phi_zero, phi_min, phi_max
        )
    ) {
        return true;
    }

    return false;
}


bool fba::Hardware::thetaphi_to_xy(
        fbg::dpair & position,
        fbg::dpair const & center, double const & theta, double const & phi,
        double theta_arm, double phi_arm, double theta_zero, double phi_zero,
        double theta_min, double phi_min, double theta_max, double phi_max,
        bool ignore_range
    ) const {
    // Note:
    //
    // The desimodel {theta,phi}_pos values are relative to
    // {theta,phi}_offset.
    //
    // However, the thetaphi_to_xy() routine (as well as other
    // routines that take theta,phi values) here expects "absolute"
    // theta,phi values -- theta angles in the x,y focal-plane frame;
    // we therefore need to apply theta,phi_offset values if we're
    // passing in theta,phi values from desimodel.
    //
    // This is way more important for theta, because desimodel has
    // applied the petal rotations to the theta_offsets -- so it
    // includes both rotation of the device and the petal position --
    // the theta_offsets span the full range of [-pi, +pi].
    // The phi_offset values are all clustered around zero.

    // Check that requested angles are in range.
    if (
        _outside_theta_phi_range(
            theta, phi, theta_zero, theta_min, theta_max, phi_zero, phi_min, phi_max
        )
    ) {
        if (!ignore_range)
            return true;
    }

    double x_theta = center.first + theta_arm * ::cos(theta);
    double y_theta = center.second + theta_arm * ::sin(theta);

    position.first = x_theta + phi_arm * ::cos(theta + phi);
    position.second = y_theta + phi_arm * ::sin(theta + phi);

    return false;
}


bool fba::Hardware::move_positioner_xy(
        fbg::shape & shptheta, fbg::shape & shpphi,
        fbg::dpair const & center, fbg::dpair const & position,
        double theta_arm, double phi_arm, double theta_zero, double phi_zero,
        double theta_min, double phi_min, double theta_max, double phi_max
        ) const {

    double phi;
    double theta;
    bool fail = xy_to_thetaphi(
        theta, phi, center, position, theta_arm, phi_arm, theta_zero,
        phi_zero, theta_min, phi_min, theta_max, phi_max
    );
    if (fail) {
        return true;
    }
    return move_positioner_thetaphi(
        shptheta, shpphi, center, theta, phi,
        theta_arm, phi_arm, theta_zero, phi_zero,
        theta_min, phi_min, theta_max, phi_max
    );
}


bool fba::Hardware::position_xy_bad(int32_t loc, fbg::dpair const & xy) const {
    double phi;
    double theta;
    if ((state.at(loc) & FIBER_STATE_STUCK) || (state.at(loc) & FIBER_STATE_BROKEN)) {
        // This positioner is stuck or has a broken fiber.  We cannot position it.
        return true;
    }
    bool fail = xy_to_thetaphi(
        theta, phi,
        loc_pos_curved_mm.at(loc),
        xy,
        loc_theta_arm.at(loc),
        loc_phi_arm.at(loc),
        loc_theta_offset.at(loc),
        loc_phi_offset.at(loc),
        loc_theta_min.at(loc),
        loc_phi_min.at(loc),
        loc_theta_max.at(loc),
        loc_phi_max.at(loc)
    );
    return fail;
}


bool fba::Hardware::loc_position_xy(
    int32_t loc, fbg::dpair const & xy, fbg::shape & shptheta,
    fbg::shape & shpphi) const {

    if ((state.at(loc) & FIBER_STATE_STUCK) || (state.at(loc) & FIBER_STATE_BROKEN)) {
        // This positioner is stuck or has a broken fiber.  We cannot move it to a
        // different X/Y location.
        return true;
    }

    // Start from exclusion polygon for this location.
    shptheta = loc_theta_excl.at(loc);
    shpphi = loc_phi_excl.at(loc);

    // std::cout << "loc_position_xy start" << std::endl;
    // shptheta.print();
    // shpphi.print();

    bool failed = move_positioner_xy(
        shptheta, shpphi,
        loc_pos_curved_mm.at(loc),
        xy,
        loc_theta_arm.at(loc),
        loc_phi_arm.at(loc),
        loc_theta_offset.at(loc),
        loc_phi_offset.at(loc),
        loc_theta_min.at(loc),
        loc_phi_min.at(loc),
        loc_theta_max.at(loc),
        loc_phi_max.at(loc)
    );

    return failed;
}


bool fba::Hardware::loc_position_thetaphi(
    int32_t loc, double theta, double phi, fbg::shape & shptheta,
    fbg::shape & shpphi, bool ignore_range) const {

    // Start from exclusion polygon for this location.
    shptheta = loc_theta_excl.at(loc);
    shpphi = loc_phi_excl.at(loc);

    bool failed = move_positioner_thetaphi(
        shptheta, shpphi,
        loc_pos_curved_mm.at(loc),
        theta, phi,
        loc_theta_arm.at(loc),
        loc_phi_arm.at(loc),
        loc_theta_offset.at(loc),
        loc_phi_offset.at(loc),
        loc_theta_min.at(loc),
        loc_phi_min.at(loc),
        loc_theta_max.at(loc),
        loc_phi_max.at(loc),
        ignore_range);

    return failed;
}


bool fba::Hardware::collide_xy(int32_t loc1, fbg::dpair const & xy1,
                               int32_t loc2, fbg::dpair const & xy2) const {

    fbg::shape shptheta1(loc_theta_excl.at(loc1));
    fbg::shape shpphi1(loc_phi_excl.at(loc1));
    bool failed1 = loc_position_xy(loc1, xy1, shptheta1, shpphi1);
    if (failed1) {
        // A positioner movement failure means that the angles needed to reach
        // the X/Y position are out of range.  While not strictly a collision,
        // it still means that we can't accept this positioner configuration.
        return true;
    }

    fbg::shape shptheta2(loc_theta_excl.at(loc2));
    fbg::shape shpphi2(loc_phi_excl.at(loc2));
    bool failed2 = loc_position_xy(loc2, xy2, shptheta2, shpphi2);
    if (failed2) {
        // A positioner movement failure means that the angles needed to reach
        // the X/Y position are out of range.  While not strictly a collision,
        // it still means that we can't accept this positioner configuration.
        return true;
    }

    // We were able to move positioners into place.  Now check for
    // intersections.

    if (fbg::intersect(shpphi1, shpphi2)) {
        return true;
    }
    if (fbg::intersect(shptheta1, shpphi2)) {
        return true;
    }
    if (fbg::intersect(shptheta2, shpphi1)) {
        return true;
    }

    return false;
}


bool fba::Hardware::collide_xy_edges(
        int32_t loc, fbg::dpair const & xy
    ) const {

    fbg::shape shptheta(loc_theta_excl.at(loc));
    fbg::shape shpphi(loc_phi_excl.at(loc));
    bool failed = loc_position_xy(loc, xy, shptheta, shpphi);
    if (failed) {
        // A positioner movement failure means that the angles needed to reach
        // the X/Y position are out of range.  While not strictly a collision,
        // it still means that we can't accept this positioner configuration.
        return true;
    }

    // We were able to move positioner into place.  Now check for
    // intersections with the GFA and petal boundaries.

    fbg::shape shpgfa(loc_gfa_excl.at(loc));
    fbg::shape shppetal(loc_petal_excl.at(loc));

    // The central body (theta arm) should never hit the GFA or petal edge,
    // so we only need to check the phi arm.

    if (fbg::intersect(shpphi, shpgfa)) {
        return true;
    }
    if (fbg::intersect(shpphi, shppetal)) {
        return true;
    }

    return false;
}


bool fba::Hardware::collide_thetaphi(
        int32_t loc1, double theta1, double phi1,
        int32_t loc2, double theta2, double phi2) const {

    fbg::shape shptheta1(loc_theta_excl.at(loc1));
    fbg::shape shpphi1(loc_phi_excl.at(loc1));
    bool failed1 = loc_position_thetaphi(loc1, theta1, phi1, shptheta1,
                                         shpphi1);
    if (failed1) {
        // A positioner movement failure means that the angles needed to reach
        // the X/Y position are out of range.  While not strictly a collision,
        // it still means that we can't accept this positioner configuration.
        return true;
    }

    fbg::shape shptheta2(loc_theta_excl.at(loc2));
    fbg::shape shpphi2(loc_phi_excl.at(loc2));
    bool failed2 = loc_position_thetaphi(loc2, theta2, phi2, shptheta2,
                                         shpphi2);
    if (failed2) {
        // A positioner movement failure means that the angles needed to reach
        // the X/Y position are out of range.  While not strictly a collision,
        // it still means that we can't accept this positioner configuration.
        return true;
    }

    // We were able to move positioners into place.  Now check for
    // intersections.

    if (fbg::intersect(shpphi1, shpphi2)) {
        return true;
    }
    if (fbg::intersect(shptheta1, shpphi2)) {
        return true;
    }
    if (fbg::intersect(shptheta2, shpphi1)) {
        return true;
    }

    return false;
}


bool fba::Hardware::collide_xy_thetaphi(
        int32_t loc1, fbg::dpair const & xy1,
        int32_t loc2, double theta2, double phi2,
        bool ignore_thetaphi_range) const {

    fbg::shape shptheta1(loc_theta_excl.at(loc1));
    fbg::shape shpphi1(loc_phi_excl.at(loc1));
    bool failed1 = loc_position_xy(loc1, xy1, shptheta1, shpphi1);
    if (failed1) {
        return true;
    }

    fbg::shape shptheta2(loc_theta_excl.at(loc2));
    fbg::shape shpphi2(loc_phi_excl.at(loc2));
    bool failed2 = loc_position_thetaphi(loc2, theta2, phi2, shptheta2,
                                         shpphi2, ignore_thetaphi_range);
    if (failed2 && !ignore_thetaphi_range) {
        return true;
    }

    if (fbg::intersect(shpphi1, shpphi2)) {
        return true;
    }
    if (fbg::intersect(shptheta1, shpphi2)) {
        return true;
    }
    if (fbg::intersect(shptheta2, shpphi1)) {
        return true;
    }

    return false;
}


std::vector <std::pair <bool, std::pair <fbg::shape, fbg::shape> > >
    fba::Hardware::loc_position_xy_multi(
        std::vector <int32_t> const & loc,
        std::vector <fbg::dpair> const & xy, int threads) const {
    size_t nlc = loc.size();
    std::vector <std::pair <bool, std::pair <fbg::shape, fbg::shape> > >
        result(nlc);

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
        fbg::shape & shptheta = result[f].second.first;
        fbg::shape & shpphi = result[f].second.second;
        result[f].first = loc_position_xy(loc[f], xy[f], shptheta, shpphi);
    }
    return result;
}


std::vector <std::pair <bool, std::pair <fbg::shape, fbg::shape> > >
    fba::Hardware::loc_position_thetaphi_multi(
        std::vector <int32_t> const & loc,
        std::vector <double> const & theta,
        std::vector <double> const & phi,
        int threads) const {
    size_t nlc = loc.size();
    std::vector <std::pair <bool, std::pair <fbg::shape, fbg::shape> > >
        result(nlc);

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

    #pragma omp parallel for schedule(static) default(none) shared(nlc, loc, theta, phi, result) num_threads(run_threads)
    for (size_t f = 0; f < nlc; ++f) {
        fbg::shape & shptheta = result[f].second.first;
        fbg::shape & shpphi = result[f].second.second;
        result[f].first = loc_position_thetaphi(loc[f], theta[f], phi[f],
                                                shptheta, shpphi);
    }
    return result;
}


std::vector <bool> fba::Hardware::check_collisions_xy(
    std::vector <int32_t> const & loc,
    std::vector <fbg::dpair> const & xy, int threads) const {

    size_t nlc = loc.size();

    auto fpos = loc_position_xy_multi(loc, xy, threads);

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

        bool failed1 = fpos[loc_indx.at(flow)].first;
        fbg::shape const & shptheta1 = fpos[loc_indx.at(flow)].second.first;
        fbg::shape const &shpphi1 = fpos[loc_indx.at(flow)].second.second;
        bool failed2 = fpos[loc_indx.at(fhigh)].first;
        fbg::shape const & shptheta2 = fpos[loc_indx.at(fhigh)].second.first;
        fbg::shape const & shpphi2 = fpos[loc_indx.at(fhigh)].second.second;

        bool hit = false;
        if (failed1 || failed2) {
            hit = true;
        } else if (fbg::intersect(shpphi1, shpphi2)) {
            hit = true;
        } else if (fbg::intersect(shptheta1, shpphi2)) {
            hit = true;
        } else if (fbg::intersect(shptheta2, shpphi1)) {
            hit = true;
        }
        if (hit) {
            #pragma omp critical
            {
                result[loc_indx[flow]] = true;
                result[loc_indx[fhigh]] = true;
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

    auto fpos = loc_position_thetaphi_multi(loc, theta, phi, threads);

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

    #pragma omp parallel for schedule(static) default(none) shared(npairs, checkpairs, fpos, result, loc_indx, theta, phi) num_threads(run_threads)
    for (size_t p = 0; p < npairs; ++p) {
        int32_t flow = checkpairs[p].first;
        int32_t fhigh = checkpairs[p].second;

        bool failed1 = fpos[loc_indx.at(flow)].first;
        fbg::shape const & shptheta1 = fpos[loc_indx.at(flow)].second.first;
        fbg::shape const & shpphi1 = fpos[loc_indx.at(flow)].second.second;
        bool failed2 = fpos[loc_indx.at(fhigh)].first;
        fbg::shape const & shptheta2 = fpos[loc_indx.at(fhigh)].second.first;
        fbg::shape const & shpphi2 = fpos[loc_indx.at(fhigh)].second.second;

        bool hit = false;
        if (failed1 || failed2) {
            hit = true;
        } else if (fbg::intersect(shpphi1, shpphi2)) {
            hit = true;
        } else if (fbg::intersect(shptheta1, shpphi2)) {
            hit = true;
        } else if (fbg::intersect(shptheta2, shpphi1)) {
            hit = true;
        }
        if (hit) {
            #pragma omp critical
            {
                result[loc_indx[flow]] = true;
                result[loc_indx[fhigh]] = true;
            }
        }
    }
    return result;
}
