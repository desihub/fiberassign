// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef UTILS_H
#define UTILS_H

#include <cstdint>
#include <cmath>

#include <iostream>
#include <chrono>
#include <memory>
#include <exception>
#include <vector>
#include <array>
#include <map>


namespace fiberassign {

#include <_version.h>

// Some helpful typedefs used when we want to access elements
// of a vector in some order sorted by a weight (both ascending
// and descending).

typedef std::pair <double, size_t> weight_index;

struct weight_compare {
    // Define this method here so that it is inline.
    bool operator() (weight_index const & lhs,
                     weight_index const & rhs) const {
        if (lhs.first == rhs.first) {
            return lhs.second < rhs.second;
        } else {
            return lhs.first > rhs.first;
        }
    }
};

struct weight_rcompare {
    // Define this method here so that it is inline.
    bool operator() (weight_index const & lhs,
                     weight_index const & rhs) const {
        if (lhs.first == rhs.first) {
            return lhs.second > rhs.second;
        } else {
            return lhs.first < rhs.first;
        }
    }
};

typedef std::pair <int32_t, int32_t> fiber_loc;

struct fiber_loc_compare {
    bool operator() (fiber_loc const & lhs,
                     fiber_loc const & rhs) const {
        return lhs.first < rhs.first;
    }
};

typedef std::pair <int32_t, int32_t> tile_loc;

struct tile_loc_compare {
    bool operator() (tile_loc const & lhs,
                     tile_loc const & rhs) const {
        if (lhs.first == rhs.first) {
            return lhs.second < rhs.second;
        } else {
            return lhs.first < rhs.first;
        }
    }
};

// These helpers are used by the header-only HTM tree and KD tree classes,
// and we don't want to modify that code for now.

void myexit (int const & i);

void myexception (std::exception const & e);

void myexception (std::exception & e);

// Simple environment control

class Environment {

    public :

        // Singleton access
        static Environment & get();

        int max_threads();
        void set_threads(int nthread);
        int current_threads();

    private :

        // This class is a singleton- constructor is private.
        Environment();

        int max_threads_;
        int cur_threads_;
};


// Simple class to help with timing parts of the code.

class Timer {

    typedef std::chrono::high_resolution_clock::time_point time_point;

    public :
        typedef std::shared_ptr <Timer> pshr;

        Timer();
        void start();
        void stop();
        void clear();
        double seconds() const;
        void report(char const * message);
        bool is_running() const;

    private :
        double total_;
        time_point start_;
        time_point stop_;
        bool running_;
        size_t calls_;

};


class GlobalTimers {

    public :

        // Singleton access
        static GlobalTimers & get();

        void start(std::string const & name);
        void stop(std::string const & name);
        double seconds(std::string const & name) const;
        bool is_running(std::string const & name) const;

        void stop_all();

        void report();

    private :

        // This class is a singleton- constructor is private.
        GlobalTimers();

        // The timer data
        std::map <std::string, Timer> data;
};


enum class log_level {
    none=0,    ///< Undefined
    debug=1,   ///< Debug
    info=2,    ///< Info
    warning=3, ///< Warning
    error=4,   ///< Error
    critical=5 ///< Critical
};

class Logger {

    public :

        // Singleton access
        static Logger & get();

        void debug_tfg(int32_t tile, int32_t loc, int64_t target,
            char const * msg);

        void debug(char const * msg);
        void info(char const * msg);
        void warning(char const * msg);
        void error(char const * msg);
        void critical(char const * msg);

        bool extra_debug();

    private :

        // This class is a singleton- constructor is private.
        Logger();

        log_level level_;
        std::string prefix_;
        int32_t debug_tile_;
        int32_t debug_loc_;
        int64_t debug_target_;
        bool debug_all_;
        bool extra_;

};


// This namespace is for all the geometry helper functions.
namespace geom {

typedef std::pair <double, double> dpair;

double sq(double const & A);

double sq(double const & A, double const & B);

double sq(dpair const & A);

double sq(dpair const & A, dpair const & B);

double norm(double const & A, double const & B);

double norm(dpair const & A);

double dist(dpair const & c1, dpair const & c2);

double scalar_prod(dpair p, dpair q, dpair d);

dpair cos_sin_angle(dpair P);

dpair sum_angles(dpair t, dpair a);

int orientation(dpair const & p, dpair const & q, dpair const & r);

void rot_pt(dpair & A, dpair const & ax, dpair const & angle);


class circle {

    public:
        typedef std::shared_ptr <circle> pshr;

        circle();

        circle(dpair const & cent, double const & rad);

        // Translation by t
        void transl(dpair const & t);

        // Rotation of the angle def by t (cos,sin) around axis axis
        void rotation(dpair const & t, dpair const & axis);

        void print() const;

        // Xmin, Xmax, Ymin, Ymax of the element
        void limits(std::array <double, 4> & lims) const;

        dpair center;

        double radius;

};

typedef std::vector <circle> circle_list;


class segments {

    public:

        typedef std::shared_ptr <segments> pshr;

        segments();

        segments(std::vector <dpair> const & pnts);

        // Translation by t
        void transl(dpair const & t);

        // Rotation of the angle def by t (cos,sin) around axis axis
        void rotation(dpair const & t, dpair const & axis);

        void print() const;

        // Xmin, Xmax, Ymin, Ymax of the element
        void limits(std::array <double, 4> & lims) const;

        std::vector <dpair> points;
};

typedef std::vector <segments> segments_list;


class shape {

    public:

        typedef std::shared_ptr <shape> pshr;

        shape();

        shape(shape const & other);

        shape(dpair const & ax, circle_list const & circs,
              segments_list const & segs);

        shape & operator=(const shape& other);

        ~shape();

        // Translation
        void transl(dpair const & t);

        // Rotation aroud polygon's axis
        void rotation (dpair const & t);

        // Rotation around origin
        void rotation_origin(dpair const & t);

        // Rotation about another anchor point
        void rotation_anchor(dpair const & t, dpair const & anchor);

        void print() const;

        // Xmin, Xmax, Ymin, Ymax of the polygon
        std::array <double, 4> limits() const;

        dpair axis;

        circle_list circle_data;

        segments_list segments_data;

};


// Intersection between combinations of segments and circles.

bool intersect_segment(dpair const & p1, dpair const & q1,
                       dpair const & p2, dpair const & q2);

bool intersect_circle(dpair const & c1, double const & r1,
                      dpair const & c2, double const & r2);

bool intersect_seg_circ(dpair const & A, dpair const & B,
                        dpair const & center, double const & rad);

bool intersect(shape const & A, shape const & B);


}


}
#endif
