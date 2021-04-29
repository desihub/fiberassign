// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef HARDWARE_H
#define HARDWARE_H

#include <cstdint>
#include <cmath>

#include <string>
#include <vector>
#include <memory>
#include <map>

#include <utils.h>

namespace fbg = fiberassign::geom;


namespace fiberassign {

// Fiber states.  OK is zero, any other bits are bad and
// indicate the reason.

#define FIBER_STATE_OK 0
#define FIBER_STATE_UNASSIGNED 1
#define FIBER_STATE_STUCK 2
#define FIBER_STATE_BROKEN 4
#define FIBER_STATE_RESTRICT 8


class Hardware : public std::enable_shared_from_this <Hardware> {

    public :

        typedef std::shared_ptr <Hardware> pshr;

        Hardware(
            std::string const & timestr,
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
            std::vector <fbg::shape> const & excl_petal
        );

        std::string time() const;

        double radial_ang2dist_CS5 (double const & theta_rad) const;

        double radial_dist2ang_CS5 (double const & dist_mm) const;

        double radial_ang2dist_curved (double const & theta_rad) const;

        double radial_dist2ang_curved (double const & arc_mm) const;

        fbg::dpair radec2xy(double const & tilera,
                            double const & tiledec,
                            double const & tiletheta,
                            double const & ra,
                            double const & dec,
                            bool use_CS5) const;

        void radec2xy_multi(double const & tilera, double const & tiledec,
                            double const & tiletheta,
                            std::vector <double> const & ra,
                            std::vector <double> const & dec,
                            std::vector <std::pair <double, double> >
                                & xy, bool use_CS5, int threads = 0) const;

        fbg::dpair xy2radec(double const & tilera, double const & tiledec,
                            double const & tiletheta,
                            double const & x_mm, double const & y_mm,
                            bool use_CS5) const;

        void xy2radec_multi(double const & tilera, double const & tiledec,
                            double const & tiletheta,
                            std::vector <double> const & x_mm,
                            std::vector <double> const & y_mm,
                            std::vector <std::pair <double, double> >
                                & radec, bool use_CS5, int threads = 0) const;

        bool move_positioner_thetaphi(
            fbg::shape & shptheta, fbg::shape & shpphi,
            fbg::dpair const & center, double theta, double phi,
            double theta_arm, double phi_arm, double theta_zero,
            double phi_zero, double theta_min, double phi_min,
            double theta_max, double phi_max,
            bool ignore_range=false) const;

        bool xy_to_thetaphi(
                double & theta, double & phi,
                fbg::dpair const & center, fbg::dpair const & position,
                double theta_arm, double phi_arm, double theta_zero,
                double phi_zero, double theta_min, double phi_min,
                double theta_max, double phi_max) const;

        bool thetaphi_to_xy(
            fbg::dpair & position,
            fbg::dpair const & center, double const & theta, double const & phi,
            double theta_arm, double phi_arm, double theta_zero, double phi_zero,
            double theta_min, double phi_min, double theta_max, double phi_max,
            bool ignore_range=false) const;

        bool position_xy_bad(int32_t loc, fbg::dpair const & xy) const;

        bool move_positioner_xy(
            fbg::shape & shptheta, fbg::shape & shpphi,
            fbg::dpair const & center, fbg::dpair const & position,
            double theta_arm, double phi_arm, double theta_zero,
            double phi_zero, double theta_min, double phi_min,
            double theta_max, double phi_max) const;

        bool loc_position_xy(
            int32_t loc, fbg::dpair const & xy, fbg::shape & shptheta,
            fbg::shape & shpphi) const;

        bool loc_position_thetaphi(
            int32_t loc, double theta, double phi, fbg::shape & shptheta,
            fbg::shape & shpphi, bool ignore_range=false) const;

        bool collide_xy(int32_t loc1, fbg::dpair const & xy1,
                        int32_t loc2, fbg::dpair const & xy2) const;

        bool collide_xy_edges(int32_t loc, fbg::dpair const & xy) const;

        bool collide_thetaphi(
            int32_t loc1, double theta1, double phi1,
            int32_t loc2, double theta2, double phi2) const;

        bool collide_xy_thetaphi(int32_t loc1, fbg::dpair const & xy1,
                                 int32_t loc2, double theta2, double phi2) const;

        std::vector <std::pair <bool, std::pair <fbg::shape, fbg::shape> > >
        loc_position_xy_multi(
            std::vector <int32_t> const & loc,
            std::vector <fbg::dpair> const & xy, int threads) const;

        std::vector <std::pair <bool, std::pair <fbg::shape, fbg::shape> > >
        loc_position_thetaphi_multi(
            std::vector <int32_t> const & loc,
            std::vector <double> const & theta,
            std::vector <double> const & phi,
            int threads) const;

        std::vector <bool> check_collisions_xy(
            std::vector <int32_t> const & loc,
            std::vector <fbg::dpair> const & xy, int threads) const;

        std::vector <bool> check_collisions_thetaphi(
            std::vector <int32_t> const & loc,
            std::vector <double> const & theta,
            std::vector <double> const & phi, int threads) const;

        // Get the platescale data
        std::vector <double> platescale_radius_mm() const;
        std::vector <double> platescale_theta_deg() const;

        // Get the radial arc length S(R)
        std::vector <double> radial_arclen() const;

        // Get the Locations for a particular device type
        std::vector <int32_t> device_locations(std::string const & type) const;

        // The (constant) total number of locations.
        int32_t nloc;

        // The number of petals
        int32_t npetal;

        // The number of slitblocks
        int32_t nslitblock;

        // The number of science fibers (device == "POS") per petal.
        int32_t nfiber_petal;

        // The full focalplane radius on the sky in degrees
        double focalplane_radius_deg;

        // The buffer region size to subtract from the total arm length
        // when considering available targets.
        double patrol_buffer_mm;

        // Locations
        std::vector <int32_t> locations;

        // Location to positioner centers in mm.  For the CS5 case these are X / Y
        // offsets in the tangent plane projection.
        std::map <int32_t, std::pair <double, double> > loc_pos_cs5_mm;
        std::map <int32_t, std::pair <double, double> > loc_pos_curved_mm;

        // Location to petal
        std::map <int32_t, int32_t> loc_petal;

        // petal to locations
        std::map <int32_t, std::vector <int32_t> > petal_locations;

        // Location to device
        std::map <int32_t, int32_t> loc_device;
        std::map <int32_t, std::string> loc_device_type;

        // Location to fiber ID
        std::map <int32_t, int32_t> loc_fiber;

        // Location to slitblock
        std::map <int32_t, int32_t> loc_slitblock;

        // Location to blockfiber
        std::map <int32_t, int32_t> loc_blockfiber;

        // The radius on the focalplane for considering which locations are
        // neighbors.
        double neighbor_radius_mm;

        // The current fiber state.
        std::map <int32_t, int32_t> state;

        // The neighbors of each device location.
        std::map <int32_t, std::vector <int32_t> > neighbors;

        // Whether the location should be considered for petal
        // boundary collisions
        std::map <int32_t, bool> petal_edge;

        // Whether the location should be considered for GFA
        // boundary collisions
        std::map <int32_t, bool> gfa_edge;

        // The positioner information:

        // The theta zero-points and range for each location.
        std::map <int32_t, double> loc_theta_offset;
        std::map <int32_t, double> loc_theta_min;
        std::map <int32_t, double> loc_theta_max;
        std::map <int32_t, double> loc_theta_pos;
        std::map <int32_t, double> loc_theta_arm;

        // The phi zero-points and range for each location.
        std::map <int32_t, double> loc_phi_offset;
        std::map <int32_t, double> loc_phi_min;
        std::map <int32_t, double> loc_phi_max;
        std::map <int32_t, double> loc_phi_pos;
        std::map <int32_t, double> loc_phi_arm;

        // The theta arm exclusion polygon for each location, in the default
        // (theta = 0.0) position
        std::map <int32_t, fbg::shape> loc_theta_excl;

        // The phi arm exclusion polygon for each location, in the default
        // (phi = 0.0) position
        std::map <int32_t, fbg::shape> loc_phi_excl;

        // The GFA exclusion polygon for each location
        std::map <int32_t, fbg::shape> loc_gfa_excl;

        // The Petal exclusion polygon for each location
        std::map <int32_t, fbg::shape> loc_petal_excl;

    private :

        std::string timestr_;

        std::vector <double> ps_radius_;
        std::vector <double> ps_theta_;
        size_t ps_size_;

        std::vector <double> arclen_;

        bool move_positioner(fbg::shape & shptheta, fbg::shape & shpphi,
                             fbg::dpair const & center,
                             fbg::dpair const & position,
                             double theta_arm, double phi_arm,
                             double theta_zero, double phi_zero,
                             double theta_min, double phi_min,
                             double theta_max, double phi_max) const;

};

}
#endif
