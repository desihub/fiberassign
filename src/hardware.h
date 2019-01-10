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


class Hardware : public std::enable_shared_from_this <Hardware> {

    public :

        typedef std::shared_ptr <Hardware> pshr;

        Hardware(
            std::vector <int32_t> const & fiber,
            std::vector <int32_t> const & petal,
            std::vector <int32_t> const & spectro,
            std::vector <int32_t> const & location,
            std::vector <int32_t> const & slit,
            std::vector <int32_t> const & slitblock,
            std::vector <int32_t> const & blockfiber,
            std::vector <int32_t> const & device,
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

        std::pair <fbg::shape, fbg::shape> fiber_position(
            int32_t fiber_id, fbg::dpair const & xy) const;

        std::vector <std::pair <fbg::shape, fbg::shape> > fiber_position_multi(
            std::vector <int32_t> const & fiber_id,
            std::vector <fbg::dpair> const & xy, int threads = 0) const;

        std::vector <bool> check_collisions_xy(
            std::vector <int32_t> const & fiber_id,
            std::vector <fbg::dpair> const & xy, int threads = 0) const;

        std::vector <bool> check_collisions_thetaphi(
            std::vector <int32_t> const & fiber_id,
            std::vector <double> const & theta,
            std::vector <double> const & phi, int threads = 0) const;

        // The (constant) total number of fibers.
        int32_t nfiber;

        // The number of petals (spectraographs)
        int32_t npetal;

        // The number of fibers per petal.
        int32_t nfiber_petal;

        // The full focalplane radius on the sky in degrees
        double focalplane_radius_deg;

        // fiber IDs
        std::vector <int32_t> fiber_id;

        // fiber ID to positioner centers in mm
        std::map <int32_t, std::pair <double, double> > fiber_pos_xy_mm;
        std::map <int32_t, double> fiber_pos_z_mm;

        // fiber ID to petal
        std::map <int32_t, int32_t> fiber_petal;

        // petal to fiber IDs
        std::map <int32_t, std::vector <int32_t> > petal_fibers;

        // fiber ID to device
        std::map <int32_t, int32_t> fiber_device;

        // fiber ID to location
        std::map <int32_t, int32_t> fiber_location;

        // fiber ID to spectrograph
        std::map <int32_t, int32_t> fiber_spectro;

        // fiber ID to positioner q/s
        std::map <int32_t, double> fiber_pos_q_deg;
        std::map <int32_t, double> fiber_pos_s_mm;

        // fiber ID to slit
        std::map <int32_t, int32_t> fiber_slit;

        // fiber ID to slitblock
        std::map <int32_t, int32_t> fiber_slitblock;

        // fiber ID to blockfiber
        std::map <int32_t, int32_t> fiber_blockfiber;

        // patrol radius in mm
        double patrol_mm;

        // The guaranteed collision distance for fibers in mm
        double collide_mm;

        // The average collision distance for fibers in mm
        double collide_avg_mm;

        // The guaranteed collision avoidance for fibers in mm
        double no_collide_mm;

        // The radius on the focalplane for considering which fibers are
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

        // The neighbors of each fiber positioner
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


//
// // Can contain information on geometry of fiber one day (hardcoded for now) not
// // used yet
// class PosP {
//     public:
//         double r1, r2;
//         PosP (double r1, double r2);
// };
//
// // Rotation of A (angle t) aroud axis axis
// void rot_pt (dpair & A, const dpair & ax, const dpair & angle);
//
// // Intersection between segments p1-q1 and p2-q2
// bool intersect (const dpair & p1, const dpair & q1, const dpair & p2,
//                 const dpair & q2);
//
// // List segments, or circles
// class element {
//     public:
//         // True if segments, false if circle
//         bool is_seg;
//         // Segments
//         Dplist segs;
//         // Circle center
//         dpair O;
//         // Circle radius
//         double rad;
//         // For python plot
//         char color;
//         // For plotting
//         double transparency;
//         // For plotting
//         double radplot;
//         element ();
//         // b : is_seg
//         element (bool b);
//         // Creates circle
//         element (const dpair & O, const double & rad);
//         // Creates segment list with only the segment AB
//         element (const dpair & A, const dpair & B);
//         // Point with color
//         element (const dpair & A, char c, double transp, double rad0);
//         // Add point to segments list
//         void add (const double & a, const double & b);
//         // Add point to segments list
//         void add (const dpair & p);
//         // Translation by t
//         void transl (const dpair & t);
//         // Rotation of the angle def by t (cos,sin) around axis axis
//         void rotation (const dpair & t, const dpair & axis);
//         void print () const;
//         // Xmin, Xmax, Ymin, Ymax of the element
//         void limits (Dlist & lims) const;
// };
//
// class Elements : public std::vector <element> {};
//
// // Intersection between 2 elements
// bool intersect (const element & e1, const element & e2);
//
// class polygon {
//     public:
//         // Segments and circles
//         Elements elmts;
//         dpair axis;
//         void add (const element & el);
//         // Add all elements of p
//         void add (const polygon & p);
//         // Translation
//         void transl (const dpair & t);
//         // Rotation aroud polygon's axis
//         void rotation (const dpair & t);
//         // Rotation around O
//         void rotation_origin (const dpair & t);
//         void print () const;
//         // Change color of the polygon
//         void set_color (char c);
//         // Xmin, Xmax, Ymin, Ymax of the polygon
//         Dlist limits () const;
// };
//
// // fh = ferrule holder
// polygon create_fh ();
// // cb = central body
// polygon create_cb ();
// // Check collision
// bool collision (const polygon & p1, const polygon & p2);
// // cos sin theta and phi for a galaxy which is in A (ref to the origin)
// Dlist angles (const dpair & A, const PosP & posp);
// // G loc of galaxy, O origin. Repositions fh and cb in the good position
// // according to O, G
// void repos_cb_fh (polygon & cb, polygon & fh, const dpair & O, const dpair & G,
//                   const PosP & posp);
//
