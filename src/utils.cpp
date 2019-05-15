// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <utils.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sstream>

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fba = fiberassign;

namespace fbg = fiberassign::geom;


void fba::myexit (int const & flag) {
    std::cout.flush();
    std::cerr.flush();
    exit(flag);
    return;
}


void fba::myexception (std::exception const & e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout.flush();
    std::cerr.flush();
    exit(1);
    return;
}


void fba::myexception (std::exception & e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout.flush();
    std::cerr.flush();
    exit(1);
    return;
}


fba::Timer::Timer() {
    clear();
}


void fba::Timer::start() {
    if (! running_) {
        start_ = std::chrono::high_resolution_clock::now();
        running_ = true;
        calls_++;
    }
    return;
}


void fba::Timer::stop() {
    if (running_) {
        stop_ = std::chrono::high_resolution_clock::now();
        std::chrono::duration <double> elapsed =
            std::chrono::duration_cast <std::chrono::duration <double> >
                (stop_ - start_);
        total_ += elapsed.count();
        running_ = false;
    }
    return;
}


void fba::Timer::clear() {
    start_ = time_point();
    stop_ = time_point();
    running_ = false;
    calls_ = 0;
    total_ = 0.0;
    return;
}


double fba::Timer::seconds() const {
    if (running_) {
        throw std::runtime_error("Timer is still running!");
    }
    return total_;
}


bool fba::Timer::is_running() const {
    return running_;
}


void fba::Timer::report(char const * message) {
    double t = seconds();
    fba::Logger & logger = fba::Logger::get();
    std::ostringstream msg;
    msg.precision(2);
    msg << std::fixed << message << ":  " << t << " seconds ("
        << calls_ << " calls)";
    logger.info(msg.str().c_str());
    return;
}


fba::Environment::Environment() {
    max_threads_ = 1;
    #ifdef _OPENMP
    max_threads_ = omp_get_max_threads();
    #endif
    cur_threads_ = max_threads_;
}

fba::Environment & fba::Environment::get() {
    static fba::Environment instance;
    return instance;
}

int fba::Environment::max_threads() {
    return max_threads_;
}

int fba::Environment::current_threads() {
    return cur_threads_;
}

void fba::Environment::set_threads(int nthread) {
    if (nthread > max_threads_) {
        auto & log = fba::Logger::get();
        std::ostringstream o;
        o << "Requested number of threads (" << nthread
            << ") is greater than the maximum (" << max_threads_
            << ") using " << max_threads_ << " instead";
        log.warning(o.str().c_str());
        nthread = max_threads_;
    }
    #ifdef _OPENMP
    omp_set_num_threads(nthread);
    #endif
    cur_threads_ = nthread;
    return;
}


fba::GlobalTimers::GlobalTimers() {
    data.clear();
}


fba::GlobalTimers & fba::GlobalTimers::get() {
    static fba::GlobalTimers instance;
    return instance;
}


void fba::GlobalTimers::start(std::string const & name) {
    if (data.count(name) == 0) {
        data[name].clear();
    }
    data.at(name).start();
    return;
}


void fba::GlobalTimers::stop(std::string const & name) {
    if (data.count(name) == 0) {
        std::ostringstream o;
        o << "Cannot stop timer " << name << " which does not exist";
        throw std::runtime_error(o.str().c_str());
    }
    data.at(name).stop();
    return;
}


double fba::GlobalTimers::seconds(std::string const & name) const {
    if (data.count(name) == 0) {
        std::ostringstream o;
        o << "Cannot get seconds for timer " << name
            << " which does not exist";
        throw std::runtime_error(o.str().c_str());
    }
    return data.at(name).seconds();
}


bool fba::GlobalTimers::is_running(std::string const & name) const {
    if (data.count(name) == 0) {
        return false;
    }
    return data.at(name).is_running();
}


void fba::GlobalTimers::stop_all() {
    for (auto & tm : data) {
        tm.second.stop();
    }
    return;
}


void fba::GlobalTimers::report() {
    stop_all();
    std::vector <std::string> names;
    for (auto const & tm : data) {
        names.push_back(tm.first);
    }
    std::stable_sort(names.begin(), names.end());
    std::ostringstream msg;
    for (auto const & nm : names) {
        msg.str("");
        msg << "Global timer: " << nm;
        data.at(nm).report(msg.str().c_str());
    }
    return;
}


fba::Logger::Logger() {
    // Prefix for messages
    prefix_ = std::string("");

    // Check DESI log level:
    level_ = log_level::none;
    char * val = ::getenv("DESI_LOGLEVEL");
    if (val == NULL) {
        level_ = log_level::info;
    } else if (strncmp(val, "DEBUG", 5) == 0) {
        level_ = log_level::debug;
    } else if (strncmp(val, "INFO", 4) == 0) {
        level_ = log_level::info;
    } else if (strncmp(val, "WARNING", 7) == 0) {
        level_ = log_level::warning;
    } else if (strncmp(val, "ERROR", 5) == 0) {
        level_ = log_level::error;
    } else if (strncmp(val, "CRITICAL", 8) == 0) {
        level_ = log_level::critical;
    }

    // Check if we are further debugging the assignment process for specific
    // tile, loc, and target values.  If these environment variables are set,
    // override the log level and set it to debug.

    extra_ = false;

    debug_all_ = false;
    val = ::getenv("DESI_DEBUG_ALL");
    if (val != NULL) {
        debug_all_ = true;
        extra_ = true;
    }

    debug_tile_ = -1;
    val = ::getenv("DESI_DEBUG_TILE");
    if (val != NULL) {
        debug_tile_ = ::atoi(val);
        extra_ = true;
    }

    debug_loc_ = -1;
    val = ::getenv("DESI_DEBUG_FIBER");
    if (val != NULL) {
        debug_loc_ = ::atoi(val);
        extra_ = true;
    }

    debug_target_ = -1;
    val = ::getenv("DESI_DEBUG_TARGET");
    if (val != NULL) {
        debug_target_ = ::atol(val);
        extra_ = true;
    }

    if (extra_) {
        if (level_ > log_level::debug) {
            fprintf(stdout, "%s: environment contains extra debug options.  Forcing DESI_LOGLEVEL=DEBUG\n", prefix_.c_str());
            fflush(stdout);
            level_ = log_level::debug;
        }
    }
}


fba::Logger & fba::Logger::get() {
    static fba::Logger instance;
    return instance;
}


void fba::Logger::debug_tfg(int32_t tile, int32_t loc, int64_t target,
    char const * msg) {
    if (level_ <= log_level::debug) {
        bool doprint = false;
        if (debug_all_) {
            doprint = true;
        } else {
            if ((debug_tile_ >= 0) && (debug_tile_ == tile)) {
                doprint = true;
            }
            if ((debug_loc_ >= 0) && (debug_loc_ == loc)) {
                doprint = true;
            }
            if ((debug_target_ >= 0) && (debug_target_ == target)) {
                doprint = true;
            }
        }
        if (doprint) {
            fprintf(stdout, "%sDEBUG: %s\n", prefix_.c_str(), msg);
            fflush(stdout);
        }
    }
    return;
}

void fba::Logger::debug(char const * msg) {
    if (level_ <= log_level::debug) {
        fprintf(stdout, "%sDEBUG: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void fba::Logger::info(char const * msg) {
    if (level_ <= log_level::info) {
        fprintf(stdout, "%sINFO: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void fba::Logger::warning(char const * msg) {
    if (level_ <= log_level::warning) {
        fprintf(stdout, "%sWARNING: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void fba::Logger::error(char const * msg) {
    if (level_ <= log_level::error) {
        fprintf(stdout, "%sERROR: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void fba::Logger::critical(char const * msg) {
    if (level_ <= log_level::critical) {
        fprintf(stdout, "%sCRITICAL: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

bool fba::Logger::extra_debug() {
    return extra_;
}


// Geometry tools

double fbg::sq(double const & A) {
    return A * A;
}

double fbg::sq(double const & A, double const & B) {
    return (A * A + B * B);
}

double fbg::sq(fbg::dpair const & A) {
    return (A.first * A.first + A.second * A.second);
}

double fbg::sq(fbg::dpair const & A, fbg::dpair const & B) {
    return (fbg::sq(A.first - B.first) + fbg::sq(A.second - B.second));
}

double fbg::norm(double const & A, double const & B) {
    return ::sqrt(A * A + B * B);
}

double fbg::norm(fbg::dpair const & A) {
    return ::sqrt(A.first * A.first + A.second * A.second);
}

double fbg::dist(fbg::dpair const & c1, fbg::dpair const & c2) {
    return fbg::norm(c1.first - c2.first, c1.second - c2.second);
}

double fbg::scalar_prod(fbg::dpair p, fbg::dpair q, fbg::dpair d) {
    return ((q.first - p.first) * (d.first - p.first)
            + (q.second - p.second) * (d.second - p.second));
}

fbg::dpair fbg::cos_sin_angle(fbg::dpair P) {
    double r_inv = 1.0 / fbg::norm(P);
    return std::make_pair(P.first * r_inv, P.second * r_inv);
}

fbg::dpair fbg::sum_angles(fbg::dpair t, fbg::dpair a) {
    return std::make_pair(t.first * a.first - t.second * a.second,
                          t.second * a.first + t.first * a.second);
}


int fbg::orientation(fbg::dpair const & p, fbg::dpair const & q,
                     fbg::dpair const & r) {
    // See 10th slides from following link for derivation of the formula
    // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    double val = (q.second - p.second) * (r.first - q.first)
                 - (q.first - p.first) * (r.second - q.second);
    // clock or counterclock wise
    return (val > 0) ? 1 : 2;
}

// Rotation of A (angle t) around axis.
void fbg::rot_pt(fbg::dpair & A, fbg::dpair const & ax,
                 fbg::dpair const & angle) {
    fbg::dpair P = std::make_pair(A.first - ax.first, A.second - ax.second);
    if ((P.first == 0) && (P.second == 0)) {
        return;
    }
    double r = fbg::norm(P);
    dpair cos_sin = fbg::cos_sin_angle(P);
    dpair sum = fbg::sum_angles(cos_sin, angle);
    A.first = ax.first + r * sum.first;
    A.second = ax.second + r * sum.second;
    return;
}


fbg::circle::circle() {
    center = std::make_pair(0.0, 0.0);
    radius = 0.0;
}


fbg::circle::circle(fbg::dpair const & cent, double const & rad) {
    center = cent;
    radius = rad;
}


void fbg::circle::transl(fbg::dpair const & t) {
    center.first += t.first;
    center.second += t.second;
    return;
}


void fbg::circle::rotation(fbg::dpair const & t, fbg::dpair const & axis) {
    fbg::rot_pt(center, axis, t);
    return;
}


void fbg::circle::print() const {
    std::cout << "  circle at (" << center.first << ", " << center.second
        << ") with radius " << radius << std::endl;
    return;
}


void fbg::circle::limits(std::array <double, 4> & lims) const {
    double xmdist = center.first - radius;
    double xpdist = center.first + radius;
    double ymdist = center.second - radius;
    double ypdist = center.second + radius;
    if (xmdist < lims[0]) {
        lims[0] = xmdist;
    }
    if (ymdist < lims[2]) {
        lims[2] = ymdist;
    }
    if (xpdist > lims[1]) {
        lims[1] = xpdist;
    }
    if (ypdist > lims[3]) {
        lims[3] = ypdist;
    }
    return;
}


fbg::segments::segments() {
    points.clear();
}


fbg::segments::segments(std::vector <fbg::dpair> const & pnts) {
    points = pnts;
}


void fbg::segments::transl(fbg::dpair const & t) {
    for (auto & s : points) {
        s.first += t.first;
        s.second += t.second;
    }
    return;
}


void fbg::segments::rotation(fbg::dpair const & t, fbg::dpair const & axis) {
    for (auto & s : points) {
        fbg::rot_pt(s, axis, t);
    }
    return;
}


void fbg::segments::print() const {
    std::cout << "  segment points: ";
    for (auto & s : points) {
        std::cout << "(" << s.first << ", " << s.second
            << ") ";
    }
    std::cout << std::endl;
    return;
}


void fbg::segments::limits(std::array <double, 4> & lims) const {
    for (auto & s : points) {
        if (s.first < lims[0]) {
            lims[0] = s.first;
        }
        if (s.first > lims[1]) {
            lims[1] = s.first;
        }
        if (s.second < lims[2]) {
            lims[2] = s.second;
        }
        if (s.second > lims[3]) {
            lims[3] = s.second;
        }
    }
    return;
}


fbg::shape::shape() {
    axis = std::make_pair(0.0, 0.0);
}


fbg::shape::shape(fbg::dpair const & ax, fbg::circle_list const & circs,
                  fbg::segments_list const & segs) {
    axis = ax;
    circle_data = circs;
    segments_data = segs;
}


void fbg::shape::transl(fbg::dpair const & t) {
    for (auto & e : circle_data) {
        e.transl(t);
    }
    for (auto & e : segments_data) {
        e.transl(t);
    }
    axis.first += t.first;
    axis.second += t.second;
    return;
}


void fbg::shape::rotation (fbg::dpair const & t) {
    for (auto & e : circle_data) {
        e.rotation(t, axis);
    }
    for (auto & e : segments_data) {
        e.rotation(t, axis);
    }
    return;
}


void fbg::shape::rotation_origin(fbg::dpair const & t) {
    fbg::dpair origin = std::make_pair(0.0, 0.0);
    for (auto & e : circle_data) {
        e.rotation(t, origin);
    }
    for (auto & e : segments_data) {
        e.rotation(t, origin);
    }
    fbg::rot_pt(axis, origin, t);
    return;
}


void fbg::shape::print() const {
    std::cout << "shape:" << std::endl;
    for (auto & e : circle_data) {
        e.print();
    }
    for (auto & e : segments_data) {
        e.print();
    }
    std::cout << "  axis = " << axis.first << ", " << axis.second << std::endl;
    return;
}


std::array <double, 4> fbg::shape::limits() const {
    std::array <double, 4> lims;
    lims[0] = 1e4;
    lims[2] = 1e4;
    lims[1] = -1e4;
    lims[3] = -1e4;
    for (auto & e : circle_data) {
        e.limits(lims);
    }
    for (auto & e : segments_data) {
        e.limits(lims);
    }
    return lims;
}


bool fbg::intersect_segment(fbg::dpair const & p1, fbg::dpair const & q1,
                            fbg::dpair const & p2, fbg::dpair const & q2) {
    int orient1 = fbg::orientation(p1, q1, p2);
    int orient2 = fbg::orientation(p1, q1, q2);
    if (orient1 == orient2) {
        return false;
    }
    int orient3 = fbg::orientation(p2, q2, p1);
    int orient4 = fbg::orientation(p2, q2, q1);
    if (orient3 == orient4) {
        return false;
    }
    return true;
}


bool fbg::intersect_circle(fbg::dpair const & c1, double const & r1,
                           fbg::dpair const & c2, double const & r2) {
    bool ret = (fbg::sq(c1, c2) < fbg::sq(r1 + r2));
    return ret;
}


// Intersection of segment and circle
bool fbg::intersect_seg_circ(fbg::dpair const & A, fbg::dpair const & B,
                             fbg::dpair const & center, double const & rad) {
    double rad_sq = fbg::sq(rad);
    double Ac_sq = fbg::sq(A, center);
    double AB_Ac = fbg::scalar_prod(A, B, center);
    if (AB_Ac <= 0.0) {
        return Ac_sq < rad_sq;
    }
    double Bc_sq = fbg::sq(B, center);
    double BA_Bc = fbg::scalar_prod(B, A, center);
    if (BA_Bc <= 0.0) {
        return Bc_sq < rad_sq;
    }
    double AB_sq = fbg::sq(A, B);
    return (Ac_sq * (1.0 - fbg::sq(AB_Ac) / (Ac_sq * AB_sq) ) < rad_sq);
}


bool fbg::intersect(fbg::shape const & A, fbg::shape const & B) {
    bool test;
    for (auto const & ac : A.circle_data) {
        for (auto const & bc : B.circle_data) {
            test = fbg::intersect_circle(ac.center, ac.radius,
                                         bc.center, bc.radius);
            if (test) {
                return true;
            }
        }
        for (auto const & bs : B.segments_data) {
            size_t bnpoints = bs.points.size();
            if (bnpoints == 0) {
                continue;
            }
            for (size_t p = 0; p < bnpoints - 1; ++p) {
                test = fbg::intersect_seg_circ(bs.points[p], bs.points[p+1],
                                               ac.center, ac.radius);
                if (test) {
                    return true;
                }
            }
        }
    }
    for (auto const & as : A.segments_data) {
        size_t anpoints = as.points.size();
        if (anpoints == 0) {
            continue;
        }
        for (auto const & bc : B.circle_data) {
            for (size_t p = 0; p < anpoints - 1; ++p) {
                test = fbg::intersect_seg_circ(as.points[p], as.points[p+1],
                                               bc.center, bc.radius);
                if (test) {
                    return true;
                }
            }
        }
        for (auto const & bs : B.segments_data) {
            size_t bnpoints = bs.points.size();
            if (bnpoints == 0) {
                continue;
            }
            for (size_t p = 0; p < anpoints - 1; ++p) {
                for (size_t q = 0; q < bnpoints - 1; ++q) {
                    test = fbg::intersect_segment(as.points[p], as.points[p+1],
                                                  bs.points[q], bs.points[q+1]);
                    if (test) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}
