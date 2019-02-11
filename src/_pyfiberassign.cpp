
#include <string>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <pybind11/stl_bind.h>

#include <hardware.h>
#include <tiles.h>
#include <targets.h>
#include <assign.h>

namespace fba = fiberassign;
namespace fbg = fiberassign::geom;
namespace py = pybind11;

using ShapeContainer = py::detail::any_container<ssize_t>;


PYBIND11_MODULE(_internal, m) {
    m.doc() = R"(
    Internal wrapper around compiled fiberassign code.

    The objects returned to python are wrapped in a shared_ptr, and wrapped
    functions are designed to take shared_ptrs as arguments.  If you are
    modifying this code, it is critical that shared_ptrs are used consistently
    everywhere.
    )";

    // Wrap the valid target types

    m.attr("TARGET_TYPE_SCIENCE") = py::int_(TARGET_TYPE_SCIENCE);
    m.attr("TARGET_TYPE_STANDARD") = py::int_(TARGET_TYPE_STANDARD);
    m.attr("TARGET_TYPE_SKY") = py::int_(TARGET_TYPE_SKY);
    m.attr("TARGET_TYPE_SAFE") = py::int_(TARGET_TYPE_SAFE);

    // Wrap the fiber states

    m.attr("FIBER_STATE_OK") = py::int_(FIBER_STATE_OK);
    m.attr("FIBER_STATE_STUCK") = py::int_(FIBER_STATE_STUCK);
    m.attr("FIBER_STATE_BROKEN") = py::int_(FIBER_STATE_BROKEN);

    py::class_ <fba::Timer, fba::Timer::pshr > (m, "Timer", R"(
        Simple timer class.

        This class is just a timer that you can start / stop / clear
        and report the results.

        )")
        .def(py::init < > ())
        .def("start", &fba::Timer::start, R"(
            Start the timer.
        )")
        .def("stop", &fba::Timer::stop, R"(
            Stop the timer.
        )")
        .def("clear", &fba::Timer::clear, R"(
            Clear the timer.
        )")
        .def("is_running", &fba::Timer::is_running, R"(
            Is the timer running?
        )")
        .def("seconds",
            [](fba::Timer const & self) {
                if (self.is_running()) {
                    return -1.0;
                } else {
                    return self.seconds();
                }
            }, R"(
            Return the elapsed seconds (if stopped) else -1.
            )"
        )
        .def("report", &fba::Timer::report, R"(
            Report results of the timer.
        )")
        .def("__repr__",
            [](fba::Timer const & self) {
                std::ostringstream o;
                o.precision(2);
                o << std::fixed;
                o << "<fiberassign.Timer ";
                if (self.is_running()) {
                    o << "(still running)";
                } else {
                    double elapsed = self.seconds();
                    o << "(stopped at " << elapsed << " seconds)";
                }
                o << ">";
                return o.str();
            }
        );


    py::class_ <fba::GlobalTimers,
        std::unique_ptr<fba::GlobalTimers, py::nodelete> > (m, "GlobalTimers",
        R"(
        Global timer registry.

        This class stores timers that can be started / stopped anywhere in
        the code to accumulate the total time for different operations.
        )")
        .def("get", [](){
            return std::unique_ptr<fba::GlobalTimers, py::nodelete>
                (&fba::GlobalTimers::get());
            }
        )
        .def("start", &fba::GlobalTimers::start)
        .def("stop", &fba::GlobalTimers::stop)
        .def("seconds", &fba::GlobalTimers::seconds)
        .def("is_running", &fba::GlobalTimers::is_running)
        .def("stop_all", &fba::GlobalTimers::stop_all)
        .def("report", &fba::GlobalTimers::report);


    py::class_ <fba::Logger, std::unique_ptr<fba::Logger, py::nodelete> > (m,
        "Logger", R"(
        Simple Logging class.

        This class mimics the python logger in C++ and respects DESI_LOGLEVEL.
        )")
        .def("get", [](){
            return std::unique_ptr<fba::Logger, py::nodelete>
                (&fba::Logger::get());
            }
        )
        .def("debug", &fba::Logger::debug, R"(
            Print a DEBUG level message.
        )")
        .def("info", &fba::Logger::info, R"(
            Print an INFO level message.
        )")
        .def("warning", &fba::Logger::warning, R"(
            Print a WARNING level message.
        )")
        .def("error", &fba::Logger::error, R"(
            Print an ERROR level message.
        )")
        .def("critical", &fba::Logger::critical, R"(
            Print a CRITICAL level message.
        )");


    py::class_ <fbg::circle, fbg::circle::pshr > (m, "Circle", R"(
        A Circle.
        )")
        .def(py::init <> ())
        .def(py::init <fbg::dpair const &, double const &> ())
        .def("transl", &fbg::circle::transl, R"(
            Translate the circle.
        )")
        .def("rotation", &fbg::circle::rotation, R"(
            Apply a rotation.
        )")
        .def_readonly("center", &fbg::circle::center)
        .def_readonly("radius", &fbg::circle::radius)
        .def("__repr__",
            [](fbg::circle const & self) {
                std::ostringstream o;
                o << "<fiberassign.Circle center=(" << self.center.first
                    << "," << self.center.second << ") radius="
                    << self.radius << ">";
                return o.str();
            }
        );


    py::class_ <fbg::segments, fbg::segments::pshr > (m, "Segments", R"(
        A collection of line segments.
        )")
        .def(py::init <> ())
        .def(py::init <std::vector <fbg::dpair> const &> ())
        .def("transl", &fbg::segments::transl, R"(
            Translate the segments.
        )")
        .def("rotation", &fbg::segments::rotation, R"(
            Apply a rotation.
        )")
        .def_readonly("points", &fbg::segments::points)
        .def("__repr__",
            [](fbg::segments const & self) {
                std::ostringstream o;
                o << "<fiberassign.Segments " << self.points.size()
                    << " points>";
                return o.str();
            }
        );


    py::class_ <fbg::shape, fbg::shape::pshr > (m, "Shape", R"(
        A shape made up of circles and segments.
        )")
        .def(py::init <> ())
        .def(py::init <fbg::dpair const &, fbg::circle_list const &,
             fbg::segments_list const &> ())
        .def("transl", &fbg::shape::transl, R"(
            Translate the shape.
        )")
        .def("rotation", &fbg::shape::rotation, R"(
            Apply a rotation.
        )")
        .def("rotation_origin", &fbg::shape::rotation_origin, R"(
            Apply a rotation about the origin.
        )")
        .def_readonly("circles", &fbg::shape::circle_data)
        .def_readonly("segments", &fbg::shape::segments_data)
        .def("__repr__",
            [](fbg::shape const & self) {
                std::ostringstream o;
                o << "<fiberassign.Shape " << self.circle_data.size()
                    << " circles and " << self.segments_data.size()
                    << " segment collections>";
                return o.str();
            }
        );


    py::class_ <fba::Hardware, fba::Hardware::pshr > (m, "Hardware", R"(
        Class representing the hardware configuration of the telescope.

        Args:
            fiber (array):  array of fiber indices (int32).
            petal (array):  array of start times (int64) for each entry.
            spectro (array):  array of stop times (int64) for each entry.
            x (array):  array of fiber X coordinate centers.
            y (array):  array of fiber Y coordinate centers.
            z (array):  array of fiber Z coordinate centers.
            status (array):  array of integers containing the fiber status.

        )")
        .def(py::init <
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <std::string> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <int32_t> const & > ())
        .def_readonly("nfiber", &fba::Hardware::nfiber)
        .def_readonly("npetal", &fba::Hardware::npetal)
        .def_readonly("nfiber_petal", &fba::Hardware::nfiber_petal)
        .def_readonly("focalplane_radius_deg",
                      &fba::Hardware::focalplane_radius_deg)
        .def_readonly("fiber_id", &fba::Hardware::fiber_id)
        .def_readonly("fiber_petal", &fba::Hardware::fiber_petal)
        .def_readonly("petal_fibers", &fba::Hardware::petal_fibers)
        .def_readonly("fiber_pos_xy_mm", &fba::Hardware::fiber_pos_xy_mm)
        .def_readonly("fiber_pos_z_mm", &fba::Hardware::fiber_pos_z_mm)
        .def_readonly("fiber_pos_q_deg", &fba::Hardware::fiber_pos_q_deg)
        .def_readonly("fiber_pos_s_mm", &fba::Hardware::fiber_pos_s_mm)
        .def_readonly("fiber_device", &fba::Hardware::fiber_device)
        .def_readonly("fiber_device_type", &fba::Hardware::fiber_device_type)
        .def_readonly("fiber_location", &fba::Hardware::fiber_location)
        .def_readonly("fiber_spectro", &fba::Hardware::fiber_spectro)
        .def_readonly("fiber_slit", &fba::Hardware::fiber_slit)
        .def_readonly("fiber_slitblock", &fba::Hardware::fiber_slitblock)
        .def_readonly("fiber_blockfiber", &fba::Hardware::fiber_blockfiber)
        .def_readonly("patrol_mm", &fba::Hardware::patrol_mm)
        .def_readonly("collide_mm", &fba::Hardware::collide_mm)
        .def_readonly("collide_avg_mm", &fba::Hardware::collide_avg_mm)
        .def_readonly("no_collide_mm", &fba::Hardware::no_collide_mm)
        .def_readonly("neighbor_radius_mm", &fba::Hardware::neighbor_radius_mm)
        .def_readonly("positioner_range", &fba::Hardware::positioner_range)
        .def_readonly("state", &fba::Hardware::state)
        .def_readonly("neighbors", &fba::Hardware::neighbors)
        .def("device_fibers", &fba::Hardware::device_fibers)
        .def("radec2xy", &fba::Hardware::radec2xy)
        .def("radec2xy_multi", [](fba::Hardware & self, double tile_ra,
            double tile_dec, std::vector <double> const & ra,
            std::vector <double> const & dec, int threads) {
            std::vector <std::pair <double, double> > xy;
            self.radec2xy_multi(tile_ra, tile_dec, ra, dec, xy, threads);
            return xy;
        }, py::return_value_policy::take_ownership)
        .def("xy2radec", &fba::Hardware::xy2radec)
        .def("xy2radec_multi", [](fba::Hardware & self, double tile_ra,
            double tile_dec, std::vector <double> const & x_mm,
            std::vector <double> const & y_mm, int threads) {
            std::vector <std::pair <double, double> > radec;
            self.xy2radec_multi(tile_ra, tile_dec, x_mm, y_mm, radec, threads);
            return radec;
        }, py::return_value_policy::take_ownership)
        .def("pos_angles", &fba::Hardware::pos_angles)
        .def("collide", &fba::Hardware::collide)
        .def("fiber_position", &fba::Hardware::fiber_position)
        .def("check_collisions_xy", &fba::Hardware::check_collisions_xy)
        .def("check_collisions_thetaphi",
             &fba::Hardware::check_collisions_thetaphi)
        .def(py::pickle(
            [](fba::Hardware const & p) { // __getstate__
                int32_t nfiber = p.fiber_id.size();
                std::vector <int32_t> fid(nfiber);
                std::vector <int32_t> petal(nfiber);
                std::vector <int32_t> spectro(nfiber);
                std::vector <int32_t> location(nfiber);
                std::vector <int32_t> slit(nfiber);
                std::vector <int32_t> slitblock(nfiber);
                std::vector <int32_t> blockfiber(nfiber);
                std::vector <int32_t> device(nfiber);
                std::vector <std::string> device_type(nfiber);
                std::vector <double> x_mm(nfiber);
                std::vector <double> y_mm(nfiber);
                std::vector <double> z_mm(nfiber);
                std::vector <double> q_deg(nfiber);
                std::vector <double> s_mm(nfiber);
                std::vector <int32_t> status(nfiber);
                for (int32_t i = 0; i < nfiber; ++i) {
                    fid[i] = p.fiber_id.at(i);
                    petal[i] = p.fiber_petal.at(fid[i]);
                    spectro[i] = p.fiber_spectro.at(fid[i]);
                    location[i] = p.fiber_location.at(fid[i]);
                    slit[i] = p.fiber_slit.at(fid[i]);
                    slitblock[i] = p.fiber_slitblock.at(fid[i]);
                    blockfiber[i] = p.fiber_blockfiber.at(fid[i]);
                    device[i] = p.fiber_device.at(fid[i]);
                    device_type[i] = p.fiber_device_type.at(fid[i]);
                    x_mm[i] = p.fiber_pos_xy_mm.at(fid[i]).first;
                    y_mm[i] = p.fiber_pos_xy_mm.at(fid[i]).second;
                    z_mm[i] = p.fiber_pos_z_mm.at(fid[i]);
                    q_deg[i] = p.fiber_pos_q_deg.at(fid[i]);
                    s_mm[i] = p.fiber_pos_s_mm.at(fid[i]);
                    status[i] = p.state.at(fid[i]);
                }
                return py::make_tuple(fid, petal, spectro, location, slit,
                    slitblock, blockfiber, device, device_type, x_mm, y_mm, z_mm, q_deg,
                    s_mm, status);
            },
            [](py::tuple t) { // __setstate__
                return new fba::Hardware(
                    t[0].cast<std::vector<int32_t> >(),
                    t[1].cast<std::vector<int32_t> >(),
                    t[2].cast<std::vector<int32_t> >(),
                    t[3].cast<std::vector<int32_t> >(),
                    t[4].cast<std::vector<int32_t> >(),
                    t[5].cast<std::vector<int32_t> >(),
                    t[6].cast<std::vector<int32_t> >(),
                    t[7].cast<std::vector<int32_t> >(),
                    t[8].cast<std::vector<std::string> >(),
                    t[9].cast<std::vector<double> >(),
                    t[10].cast<std::vector<double> >(),
                    t[11].cast<std::vector<double> >(),
                    t[12].cast<std::vector<double> >(),
                    t[13].cast<std::vector<double> >(),
                    t[14].cast<std::vector<int32_t> >()
                );
            }
        ));


    py::class_ <fba::Target, fba::Target::pshr > (m, "Target", R"(
        Class representing a single target.
        )")
        .def(py::init <> ())
        .def_readwrite("id", &fba::Target::id)
        .def_readwrite("ra", &fba::Target::ra)
        .def_readwrite("dec", &fba::Target::dec)
        .def_readwrite("desi_target", &fba::Target::desi_target)
        .def_readwrite("bgs_target", &fba::Target::bgs_target)
        .def_readwrite("mws_target", &fba::Target::mws_target)
        .def_readwrite("obs_remain", &fba::Target::obs_remain)
        .def_readwrite("priority", &fba::Target::priority)
        .def_readwrite("subpriority", &fba::Target::subpriority)
        .def_readwrite("obscond", &fba::Target::obscond)
        .def_readwrite("type", &fba::Target::type)
        .def("is_science", &fba::Target::is_science)
        .def("is_standard", &fba::Target::is_standard)
        .def("is_sky", &fba::Target::is_sky)
        .def("is_safe", &fba::Target::is_safe)
        .def("__repr__",
            [](fba::Target const & tg) {
                std::ostringstream o;
                o << "<fiberassign.Target id=" << tg.id << " ra=" << tg.ra
                << " dec=" << tg.dec << " type=" << (int)tg.type
                << " priority=" << tg.priority << " subpriority="
                << tg.subpriority << " >";
                return o.str();
            }
        );

    // Define a numpy dtype wrapper around our internal "Target" class.
    PYBIND11_NUMPY_DTYPE(fba::Target, id, ra, dec, obs_remain, priority,
                         subpriority, obscond, type);


    py::class_ <fba::Targets, fba::Targets::pshr > (m, "Targets", R"(
        Class representing a list of targets.

        The target data is stored internally as a dictionary of Target
        instances.  After construction this container is empty, and then one
        uses the append method to add Targets.

        )")
        .def(py::init < > ())
        .def("append", &fba::Targets::append, R"(
            Append objects to the target list.

            Args:
                ids (array):  array of int64 target IDs.
                ras (array):  array of float64 target RA coordinates.
                decs (array):  array of float64 target DEC coordinates.
                obs_remain (array):  array of int32 number of remaining
                    observations.
                desi_target (array):  array of int64 DESI_TARGET values.
                bgs_target (array):  array of int64 BGS_TARGET values.
                mws_target (array):  array of int64 MWS_TARGET values.
                priority (array):  array of int32 values representing the target
                    class priority for each object.
                subpriority (array):  array of float64 values in [0.0, 1.0]
                    representing the priority within the target class.
                obscond (array):  array of int32 bitfields describing the
                    valid observing conditions for each target.
                type (array):  array of uint8 bitfields holding the types of
                    of each target (science, standard, etc).

        )")
        .def("ids", [](fba::Targets & self) {
            auto ntarg = self.data.size();
            // Create a numpy array to return
            py::array_t < int64_t > ret;
            ret.resize( {ntarg} );
            py::buffer_info info = ret.request();
            int64_t * raw = static_cast < int64_t * > (info.ptr);
            // Copy the target IDs into the numpy buffer
            size_t indx = 0;
            for (auto const & it : self.data) {
                raw[indx] = it.first;
                indx++;
            }
            // Sort the result
            std::sort(raw, &(raw[ntarg-1]));
            return ret;
        })
        .def("get", [](fba::Targets & self, int64_t id) {
            return self.data[id];
        }, py::return_value_policy::reference_internal)
        .def("__repr__",
            [](fba::Targets const & tgs) {
                std::ostringstream o;
                o << "<fiberassign.Targets with " << tgs.data.size() << " objects>";
                return o.str();
            }
        );


    py::class_ <fba::Tiles, fba::Tiles::pshr > (m, "Tiles", R"(
        Class representing a list of tiles.

        The tile data is stored internally as a vector of Tile objects.
        The constructor of this class takes vectors of tile properties
        that are easy to obtain in python when reading the tile files.

        Args:
            hw (Hardware):  the hardware description.
            ids (array):  array of int32 tile IDs.
            ras (array):  array of float64 tile RA coordinates.
            decs (array):  array of float64 tile DEC coordinates.
            obs (array):  array of int32 obsconditions bitfields.

        )")
        .def(py::init < fba::Hardware::pshr, std::vector <int32_t> const &,
            std::vector <double> const &, std::vector <double> const &,
            std::vector <int32_t> const & > ())
        .def_readonly("id", &fba::Tiles::id)
        .def_readonly("ra", &fba::Tiles::ra)
        .def_readonly("dec", &fba::Tiles::dec)
        .def_readonly("obscond", &fba::Tiles::obscond)
        .def_readonly("order", &fba::Tiles::order)
        .def("__repr__",
            [](fba::Tiles const & tls) {
                std::ostringstream o;
                o << "<fiberassign.Tiles with " << tls.id.size() << " tiles>";
                return o.str();
            }
        );


    py::class_ <fba::TargetTree, fba::TargetTree::pshr > (m, "TargetTree",
        R"(
        Class representing an HTM tree of targets

        This encapsulates a Hierarchical Triangular Mesh tree structure that
        makes it efficient search for targets within a radius of some point.

        Args:
            tgs (Targets):  A Targets object.
            min_tree_size (float):  (Optional) minimum node size of tree.

        )")
        .def(py::init < fba::Targets::pshr, double > (),
            py::arg(), py::arg("min_tree_size") = 0.01)
        .def("near", [](fba::TargetTree & self, double ra_deg,
            double dec_deg, double radius_rad) {
            std::vector <int64_t> result;
            self.near(ra_deg, dec_deg, radius_rad, result);
            return result;
        }, py::return_value_policy::take_ownership);

    py::class_ <fba::TargetsAvailable, fba::TargetsAvailable::pshr > (m,
        "TargetsAvailable", R"(
        Class representing the objects reachable by each fiber of each tile.

        This data structure makes it convenient and efficient to get the list
        of target IDs available to a given fiber for a given tile.

        Args:
            objs (Targets):  the objects.
            tiles (Tiles):  the tiles to consider.
            tree (TargetTree):  the HTM tree of objects.

        )")
        .def(py::init < fba::Targets::pshr, fba::Tiles::pshr,
            fba::TargetTree::pshr > ())
        .def("tile_data", &fba::TargetsAvailable::tile_data,
            py::return_value_policy::reference_internal);



    py::class_ <fba::FibersAvailable, fba::FibersAvailable::pshr > (m,
        "FibersAvailable", R"(
        Class representing the tile and fibers reachable by each target.

        This data structure makes it convenient and efficient to get the list
        of tile / fiber combinations available for a given target.

        Args:
            tgsavail (TargetsAvailable):  the targets available to each fiber.

        )")
        .def(py::init < fba::TargetsAvailable::pshr > ())
        .def("target_data", &fba::FibersAvailable::target_data);


    py::class_ <fba::Assignment, fba::Assignment::pshr > (m,
        "Assignment", R"(
        Class representing the current assignment of all fibers.

        This data structure stores the current fiber assignment and
        provides methods for manipulating that assignment

        Args:
            tgs (Targets):  the targets.
            tgsavail (TargetsAvailable):  the targets available to each fiber.
            favail (FibersAvailable):  the fibers available to each target.

        )")
        .def(py::init < fba::Targets::pshr, fba::TargetsAvailable::pshr,
             fba::FibersAvailable::pshr > ())
        .def("targets", &fba::Assignment::targets)
        .def("hardware", &fba::Assignment::hardware)
        .def("targets_avail", &fba::Assignment::targets_avail)
        .def("fibers_avail", &fba::Assignment::fibers_avail)
        .def("tiles_assigned", &fba::Assignment::tiles_assigned)
        .def("tile_fiber_target", &fba::Assignment::tile_fiber_target,
            py::return_value_policy::reference_internal)
        .def("assign_unused", &fba::Assignment::assign_unused,
             py::arg("tgtype")=TARGET_TYPE_SCIENCE,
             py::arg("max_per_petal")=-1,
             py::arg("pos_type")=std::string("POS"),
             py::arg("start_tile")=-1, py::arg("stop_tile")=-1,
             py::arg("max_standards_petal")=-1)
        .def("assign_force", &fba::Assignment::assign_force,
             py::arg("tgtype")=TARGET_TYPE_SCIENCE,
             py::arg("required_per_petal")=0,
             py::arg("start_tile")=-1, py::arg("stop_tile")=-1)
        .def("redistribute_science", &fba::Assignment::redistribute_science,
             py::arg("start_tile")=-1, py::arg("stop_tile")=-1,
             py::arg("max_standards_petal")=-1);


}
