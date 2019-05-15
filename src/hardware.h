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
#define FIBER_STATE_STUCK 1
#define FIBER_STATE_BROKEN 2
#define FIBER_STATE_SAFE 4


class Hardware : public std::enable_shared_from_this <Hardware> {

    public :

        typedef std::shared_ptr <Hardware> pshr;

        Hardware(
            std::vector <int32_t> const & location,
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
            std::vector <int32_t> const & status
        );

        double radial_ang2dist (double const & theta_rad) const;

        double radial_dist2ang (double const & dist_mm) const;

        fbg::dpair radec2xy(double const & tilera,
                                          double const & tiledec,
                                          double const & ra,
                                          double const & dec) const;

        void radec2xy_multi(double const & tilera, double const & tiledec,
                            std::vector <double> const & ra,
                            std::vector <double> const & dec,
                            std::vector <std::pair <double, double> >
                                & xy, int threads = 0) const;

        fbg::dpair xy2radec(double const & tilera, double const & tiledec,
                            double const & x_mm, double const & y_mm) const;

        void xy2radec_multi(double const & tilera, double const & tiledec,
                            std::vector <double> const & x_mm,
                            std::vector <double> const & y_mm,
                            std::vector <std::pair <double, double> >
                                & radec, int threads = 0) const;

        std::array <double, 4> pos_angles(fbg::dpair const & A,
            fbg::dpair const & posp) const;

        bool collide(fbg::dpair center1,
                     fbg::dpair position1,
                     fbg::dpair center2,
                     fbg::dpair position2) const;

        std::pair <fbg::shape, fbg::shape> loc_position(
            int32_t loc, fbg::dpair const & xy) const;

        std::vector <std::pair <fbg::shape, fbg::shape> > loc_position_multi(
            std::vector <int32_t> const & loc,
            std::vector <fbg::dpair> const & xy, int threads = 0) const;

        std::vector <bool> check_collisions_xy(
            std::vector <int32_t> const & loc,
            std::vector <fbg::dpair> const & xy, int threads = 0) const;

        std::vector <bool> check_collisions_thetaphi(
            std::vector <int32_t> const & loc,
            std::vector <double> const & theta,
            std::vector <double> const & phi, int threads = 0) const;

        // Get the Locations for a particular device type
        std::vector <int32_t> device_locations(std::string const & type) const;

        // The (constant) total number of locations.
        int32_t nloc;

        // The number of petals (spectraographs)
        int32_t npetal;

        // The number of fibers (device == "POS") per petal.
        int32_t nfiber_petal;

        // The full focalplane radius on the sky in degrees
        double focalplane_radius_deg;

        // Locations
        std::vector <int32_t> locations;

        // Location to positioner centers in mm
        std::map <int32_t, std::pair <double, double> > loc_pos_xy_mm;
        std::map <int32_t, double> loc_pos_z_mm;

        // Location to petal
        std::map <int32_t, int32_t> loc_petal;

        // petal to locations
        std::map <int32_t, std::vector <int32_t> > petal_locations;

        // Location to device
        std::map <int32_t, int32_t> loc_device;
        std::map <int32_t, std::string> loc_device_type;

        // Location to fiber ID
        std::map <int32_t, int32_t> loc_fiber;

        // Location to spectrograph
        std::map <int32_t, int32_t> loc_spectro;

        // Location to positioner q/s
        std::map <int32_t, double> loc_pos_q_deg;
        std::map <int32_t, double> loc_pos_s_mm;

        // Location to slit
        std::map <int32_t, int32_t> loc_slit;

        // Location to slitblock
        std::map <int32_t, int32_t> loc_slitblock;

        // Location to blockfiber
        std::map <int32_t, int32_t> loc_blockfiber;

        // patrol radius in mm
        double patrol_mm;

        // The guaranteed collision distance for positioners in mm
        double collide_mm;

        // The average collision distance for positioners in mm
        double collide_avg_mm;

        // The guaranteed collision avoidance for positioners in mm
        double no_collide_mm;

        // The radius on the focalplane for considering which locations are
        // neighbors.
        double neighbor_radius_mm;

        // The range of the positioner.
        fbg::dpair positioner_range;
        double pos_theta_max;
        double pos_phi_max;
        double pos_theta_max_deg;
        double pos_phi_max_deg;

        // The current fiber state.
        std::map <int32_t, int32_t> state;

        // The neighbors of each device location.
        std::map <int32_t, std::vector <int32_t> > neighbors;

        // The geometric shape of the ferrule holder
        fbg::shape ferrule_holder;

        // The geometric shape of the central body
        fbg::shape central_body;

    private :

        fbg::shape create_ferrule_holder() const;

        fbg::shape create_central_body() const;

        void reposition_fiber(fbg::shape & cb,
                              fbg::shape & fh,
                              fbg::dpair const & center,
                              fbg::dpair const & position,
                              fbg::dpair const & posp) const;

};

}
#endif
