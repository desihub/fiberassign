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
    m.attr("TARGET_TYPE_SUPPSKY") = py::int_(TARGET_TYPE_SUPPSKY);
    m.attr("TARGET_TYPE_SAFE") = py::int_(TARGET_TYPE_SAFE);

    // Wrap the fiber states

    m.attr("FIBER_STATE_OK") = py::int_(FIBER_STATE_OK);
    m.attr("FIBER_STATE_UNASSIGNED") = py::int_(FIBER_STATE_UNASSIGNED);
    m.attr("FIBER_STATE_STUCK") = py::int_(FIBER_STATE_STUCK);
    m.attr("FIBER_STATE_BROKEN") = py::int_(FIBER_STATE_BROKEN);
    m.attr("FIBER_STATE_SAFE") = py::int_(FIBER_STATE_SAFE);

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

            Returns:
                (bool): True if the timer is running, else False.

        )")
        .def("seconds",
            [](fba::Timer const & self) {
                if (self.is_running()) {
                    return -1.0;
                } else {
                    return self.seconds();
                }
            }, R"(
            Return the elapsed seconds.

            Returns:
                (float): The elapsed seconds (if timer is stopped) else -1.

        )")
        .def("report", &fba::Timer::report, py::arg("msg"), R"(
            Report results of the timer.

            Args:
                msg (str): The message to print before the timer value.

            Returns:
                None

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
            }, R"(
            Get the instance of global singleton class.
        )")
        .def("start", &fba::GlobalTimers::start, py::arg("name"), R"(
            Start the specified timer.

            If the named timer does not exist, it is first created before
            being started.

            Args:
                name (str): The name of the global timer.

            Returns:
                None
        )")
        .def("stop", &fba::GlobalTimers::stop, py::arg("name"), R"(
            Stop the specified timer.

            The timer must already exist.

            Args:
                name (str): The name of the global timer.

            Returns:
                None
        )")
        .def("seconds", &fba::GlobalTimers::seconds, py::arg("name"), R"(
            Get the elapsed time for a timer.

            The timer must be stopped.

            Args:
                name (str): The name of the global timer.

            Returns:
                (float): The elapsed time in seconds.
        )")
        .def("is_running", &fba::GlobalTimers::is_running, py::arg("name"), R"(
            Is the specified timer running?

            Args:
                name (str): The name of the global timer.

            Returns:
                (bool): True if the timer is running, else False.
        )")
        .def("stop_all", &fba::GlobalTimers::stop_all, R"(
            Stop all global timers.
        )")
        .def("report", &fba::GlobalTimers::report, R"(
            Report results of all global timers to STDOUT.
        )");


    py::class_ <fba::Logger, std::unique_ptr<fba::Logger, py::nodelete> > (m,
        "Logger", R"(
        Simple Logging class.

        This class mimics the python logger in C++ and respects DESI_LOGLEVEL.
        )")
        .def("get", [](){
            return std::unique_ptr<fba::Logger, py::nodelete>
                (&fba::Logger::get());
            }, R"(
            Get the instance of global singleton class.
        )")
        .def("debug", &fba::Logger::debug, py::arg("msg"), R"(
            Print a DEBUG level message.

            Args:
                msg (str): The message to print.

            Returns:
                None

        )")
        .def("info", &fba::Logger::info, py::arg("msg"), R"(
            Print an INFO level message.

            Args:
                msg (str): The message to print.

            Returns:
                None

        )")
        .def("warning", &fba::Logger::warning, py::arg("msg"), R"(
            Print a WARNING level message.

            Args:
                msg (str): The message to print.

            Returns:
                None

        )")
        .def("error", &fba::Logger::error, py::arg("msg"), R"(
            Print an ERROR level message.

            Args:
                msg (str): The message to print.

            Returns:
                None

        )")
        .def("critical", &fba::Logger::critical, py::arg("msg"), R"(
            Print a CRITICAL level message.

            Args:
                msg (str): The message to print.

            Returns:
                None

        )");


    py::class_ <fba::Environment,
        std::unique_ptr<fba::Environment, py::nodelete> > (
            m, "Environment", R"(
        Environment control.

        This class allows setting the threading concurrency of the compiled
        code.
        )")
        .def("get", [](){
            return std::unique_ptr<fba::Environment, py::nodelete>
                (&fba::Environment::get());
            }, R"(
            Get the instance of global singleton class.
        )")
        .def("max_threads", &fba::Environment::max_threads, R"(
            Return the maximum threads supported by the runtime environment.
        )")
        .def("current_threads", &fba::Environment::current_threads, R"(
            Return the current threading concurrency in use.
        )")
        .def("set_threads", &fba::Environment::set_threads,
            py::arg("nthread"), R"(
            Set the number of threads in use.

            Args:
                nthread (int): The number of threads to use.

            Returns:
                None

        )");


    py::class_ <fbg::circle, fbg::circle::pshr > (m, "Circle", R"(
        A Circle.

        This class represents a circle with a center and radius.
        This shape can be translated and rotated about an axis.

        Args:
            center (tuple): The (X, Y) center of the circle.
            radius (float): The radius of the circle.

        )")
        .def(py::init <> ())
        .def(py::init <fbg::dpair const &, double const &> (),
            py::arg("center"), py::arg("radius"))
        .def("transl", &fbg::circle::transl, py::arg("offset"), R"(
            Translate the circle.

            Args:
                offset (tuple): The (X, Y) offset to add to the center.

        )")
        .def("rotation", &fbg::circle::rotation, py::arg("angle"),
            py::arg("axis"), R"(
            Apply a rotation.

            Rotate the circle center by an angle about the given point.

            Args:
                angle (float): The angle of rotation.
                axis (tuple): The (X, Y) origin of the rotation.

        )")
        .def_readonly("center", &fbg::circle::center, R"(
            The center (X, Y) tuple of the circle.
        )")
        .def_readonly("radius", &fbg::circle::radius, R"(
            The radius of the circle.
        )")
        .def(py::pickle(
            [](fbg::circle const & p) { // __getstate__
                fbg::dpair cent = p.center;
                double rad = p.radius;
                return py::make_tuple(cent, rad);
            },
            [](py::tuple t) { // __setstate__
                return new fbg::circle(
                    t[0].cast<fbg::dpair>(),
                    t[1].cast<double>()
                );
            }
        ))
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

        This class represents a sequence of connected line segments.

        Args:
            points (list):  A list of (X, Y) tuples containing the points
                describing the segments.

        )")
        .def(py::init <> ())
        .def(py::init <std::vector <fbg::dpair> const &> (), py::arg("points"))
        .def("transl", &fbg::segments::transl, py::arg("offset"), R"(
            Translate the segments.

            Args:
                offset (tuple): The (X, Y) offset to add to all points.

        )")
        .def("rotation", &fbg::segments::rotation, py::arg("angle"),
            py::arg("axis"), R"(
            Apply a rotation.

            Rotate all points by an angle about the given origin.

            Args:
                angle (float): The angle of rotation.
                axis (tuple): The (X, Y) origin of the rotation.

        )")
        .def_readonly("points", &fbg::segments::points, R"(
            The list of points.
        )")
        .def(py::pickle(
            [](fbg::segments const & p) { // __getstate__
                std::vector <fbg::dpair> pts = p.points;
                return py::make_tuple(pts);
            },
            [](py::tuple t) { // __setstate__
                return new fbg::segments(
                    t[0].cast<std::vector <fbg::dpair> >()
                );
            }
        ))
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

        A Shape contains a list Circle objects and a list of Segments objects
        that can be rotated / translated together.  Each shape also has a
        center point independent of the individual pieces.

        Args:
            center (tuple): The (X, Y) center of the shape.
            circles (list): A list of Circle objects.
            segments (list): A list of Segments objects.

        )")
        .def(py::init <> ())
        .def(py::init <fbg::shape const &> (), py::arg("other")
        )
        .def(py::init <fbg::dpair const &, fbg::circle_list const &,
             fbg::segments_list const &> (), py::arg("center"),
             py::arg("circles"), py::arg("segments")
        )
        .def("transl", &fbg::shape::transl, py::arg("offset"), R"(
            Translate the shape.

            This adds an offset to the center and to all constituent
            circles and segments.

            Args:
                offset (tuple): The (X, Y) offset to add.

        )")
        .def("rotation", &fbg::shape::rotation, py::arg("angle"), R"(
            Apply a rotation about the center.

            Rotate all circles and segments about the center of the shape.

            Args:
                angle (float): The angle of rotation.

        )")
        .def("rotation_origin", &fbg::shape::rotation_origin,
            py::arg("angle"), R"(
            Apply a rotation about the origin.

            Rotate the entire shape about the origin.

            Args:
                angle (float): The angle of rotation.

        )")
        .def_readonly("axis", &fbg::shape::axis, R"(
            The axis (center).
        )")
        .def_readonly("circles", &fbg::shape::circle_data, R"(
            The list of constituent circles.
        )")
        .def_readonly("segments", &fbg::shape::segments_data, R"(
            The list of constituent segments.
        )")
        .def(py::pickle(
            [](fbg::shape const & p) { // __getstate__
                fbg::dpair axis = p.axis;
                fbg::circle_list cdata = p.circle_data;
                fbg::segments_list sdata = p.segments_data;
                return py::make_tuple(axis, cdata, sdata);
            },
            [](py::tuple t) { // __setstate__
                return new fbg::shape(
                    t[0].cast<fbg::dpair>(),
                    t[1].cast<fbg::circle_list>(),
                    t[2].cast<fbg::segments_list>()
                );
            }
        ))
        .def("__repr__",
            [](fbg::shape const & self) {
                std::ostringstream o;
                o << "<fiberassign.Shape at (" << self.axis.first
                    << ", " << self.axis.second << ") has "
                    << self.circle_data.size()
                    << " circles and " << self.segments_data.size()
                    << " segment collections:" << std::endl;
                for (auto const & c : self.circle_data) {
                    o << "  circle: (" << c.center.first
                    << ", " << c.center.second << ") rad = "
                    << c.radius << std::endl;
                }
                for (auto const & s : self.segments_data) {
                    o << "  segment points:" << std::endl;
                    for (auto const & p : s.points) {
                        o << "    (" << p.first << ", " << p.second
                            << ")" << std::endl;
                    }
                }
                o << ">";
                return o.str();
            }
        );


    py::class_ <fba::Hardware, fba::Hardware::pshr > (m, "Hardware", R"(
        Class representing the hardware configuration of the telescope.

        Args:
            timestr (str):  ISO 8601 format time string in UTC.
            location (array):  int32 array of location.
            petal (array):  int32 array of petal index.
            device (array):  int32 array of device number.
            slitblock (array):  int32 array of slitblock.
            blockfiber (array):  int32 array of blockfiber.
            fiber (array):  int32 array of fiber IDs.
            device_type (list):  List of strings of device types
                (e.g. POS, ETC).
            x (array):  location X coordinate centers in mm.
            y (array):  location Y coordinate centers in mm.
            status (array):  array of integers containing the fiber status.
            theta_offset (array):  The theta angle zero-points in degrees for
                each device.
            theta_min (array):  The theta angle minimum value in degrees
                relative to the offset.
            theta_max (array):  The theta angle maximum value in degrees
                relative to the offset.
            theta_arm (array):  The theta arm lengths in mm.
            phi_offset (array):  The phi angle zero-points in degrees for
                each device.
            phi_min (array):  The phi angle minimum value in degrees
                relative to the offset.
            phi_max (array):  The phi angle maximum value in degrees
                relative to the offset.
            phi_arm (array):  The phi arm lengths in mm.
            ps_radius (array):  The platescale radius vector in mm.
            ps_theta (array):  The platescale theta vector in degrees.
            arclen (array):  The radial arc length S(R) in mm.
            excl_theta (list):  The Shape object for the exclusion polygon
                of the theta arm of each device.
            excl_phi (list):  The Shape object for the exclusion polygon
                of the phi arm of each device.
            excl_gfa (list):  The Shape object for the exclusion polygon
                of the GFA for each device.
            excl_petal (list):  The Shape object for the exclusion polygon
                of the petal edge for each device.

        )")
        .def(py::init <
            std::string const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <int32_t> const &,
            std::vector <std::string> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <int32_t> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <double> const &,
            std::vector <fbg::shape> const &,
            std::vector <fbg::shape> const &,
            std::vector <fbg::shape> const &,
            std::vector <fbg::shape> const &> (), py::arg("timestr"),
            py::arg("location"),
            py::arg("petal"), py::arg("device"), py::arg("slitblock"),
            py::arg("blockfiber"), py::arg("fiber"), py::arg("device_type"),
            py::arg("x"), py::arg("y"), py::arg("status"),
            py::arg("theta_offset"), py::arg("theta_min"),
            py::arg("theta_max"), py::arg("theta_arm"),
            py::arg("phi_offset"), py::arg("phi_min"),
            py::arg("phi_max"), py::arg("phi_arm"),
            py::arg("ps_radius"), py::arg("ps_theta"), py::arg("arclen"),
            py::arg("excl_theta"), py::arg("excl_phi"),
            py::arg("excl_gfa"), py::arg("excl_petal")
        )
        .def_readonly("nloc", &fba::Hardware::nloc, R"(
            The number of device locations.
        )")
        .def_readonly("npetal", &fba::Hardware::npetal, R"(
            The number of petals.
        )")
        .def_readonly("nfiber_petal", &fba::Hardware::nfiber_petal, R"(
            The number of science positioners (device type POS) per petal.
        )")
        .def_readonly("focalplane_radius_deg",
                      &fba::Hardware::focalplane_radius_deg, R"(
                          The focalplane radius in degrees.
                      )")
        .def_readonly("locations", &fba::Hardware::locations, R"(
            Vector of locations.
        )")
        .def_readonly("loc_petal", &fba::Hardware::loc_petal, R"(
            Dictionary of the petal for each location.
        )")
        .def_readonly("petal_locations", &fba::Hardware::petal_locations, R"(
            Dictionary of the locations for each petal.
        )")
        .def_readonly("loc_pos_cs5_mm", &fba::Hardware::loc_pos_cs5_mm, R"(
            Dictionary of central (X, Y) position tuples for each location in CS5.
        )")
        .def_readonly("loc_pos_curved_mm", &fba::Hardware::loc_pos_curved_mm, R"(
            Dictionary of central position tuples for each location in curved
            focal surface coordinates.
        )")
        .def_readonly("loc_device", &fba::Hardware::loc_device, R"(
            Dictionary of device ID for each location
        )")
        .def_readonly("loc_device_type",
                      &fba::Hardware::loc_device_type, R"(
            Dictionary of device type (POS or ETC) for each location.
        )")
        .def_readonly("loc_fiber", &fba::Hardware::loc_fiber, R"(
            Dictionary of fiber values for each location.
        )")
        .def_readonly("loc_slitblock", &fba::Hardware::loc_slitblock, R"(
            Dictionary of slitblock values for each location.
        )")
        .def_readonly("loc_blockfiber", &fba::Hardware::loc_blockfiber, R"(
            Dictionary of blockfiber values for each location.
        )")
        .def_readonly("loc_theta_arm", &fba::Hardware::loc_theta_arm, R"(
            Dictionary of theta arm lengths for each location.
        )")
        .def_readonly("loc_theta_offset", &fba::Hardware::loc_theta_offset, R"(
            Dictionary of theta angle offsets for each location.
        )")
        .def_readonly("loc_theta_min", &fba::Hardware::loc_theta_min, R"(
            Dictionary of theta min range from offset for each location.
        )")
        .def_readonly("loc_theta_max", &fba::Hardware::loc_theta_max, R"(
            Dictionary of theta max range from offset for each location.
        )")
        .def_readonly("loc_phi_arm", &fba::Hardware::loc_phi_arm, R"(
            Dictionary of phi arm lengths for each location.
        )")
        .def_readonly("loc_phi_offset", &fba::Hardware::loc_phi_offset, R"(
            Dictionary of phi angle offsets for each location.
        )")
        .def_readonly("loc_phi_min", &fba::Hardware::loc_phi_min, R"(
            Dictionary of phi min range from offset for each location.
        )")
        .def_readonly("loc_phi_max", &fba::Hardware::loc_phi_max, R"(
            Dictionary of phi max range from offset for each location.
        )")
        .def_readonly("loc_theta_excl", &fba::Hardware::loc_theta_excl, R"(
            Dictionary of theta exclusion shapes for each location.
        )")
        .def_readonly("loc_phi_excl", &fba::Hardware::loc_phi_excl, R"(
            Dictionary of phi exclusion shapes for each location.
        )")
        .def_readonly("loc_gfa_excl", &fba::Hardware::loc_gfa_excl, R"(
            Dictionary of GFA exclusion shapes for each location.
        )")
        .def_readonly("loc_petal_excl", &fba::Hardware::loc_petal_excl, R"(
            Dictionary of petal exclusion shapes for each location.
        )")
        .def_readonly("neighbor_radius_mm",
                      &fba::Hardware::neighbor_radius_mm, R"(
            Radius for considering locations as neighbors.
        )")
        .def_readonly("patrol_buffer_mm",
                      &fba::Hardware::patrol_buffer_mm, R"(
            Buffer to subtract from full patrol radius when considering targets.
        )")
        .def_readonly("state", &fba::Hardware::state, R"(
            Dictionary of fiber state for each location.
        )")
        .def_readonly("neighbors", &fba::Hardware::neighbors, R"(
            Dictionary of neighbor IDs for each location.
        )")
        .def("device_locations", &fba::Hardware::device_locations,
            py::arg("type"), R"(
            Dictionary of locations for each device type (POS or ETC).
        )")
        .def("time", &fba::Hardware::time, R"(
            Return the time used when loading the focalplane model.

            Returns:
                (str): the focalplane model time.

        )")
        .def("radial_ang2dist_CS5", &fba::Hardware::radial_ang2dist_CS5,
            py::arg("theta_rad"), R"(
            Covert the angle from the origin to CS5 distance in mm.

            This uses the radial platescale to convert the angle to mm in the
            tangent plane (CS5) coordinates.

            Args:
                theta_rad (float): Theta angle in radians.

            Returns:
                (float): the distance in mm.

        )")
        .def("radial_dist2ang_CS5", &fba::Hardware::radial_dist2ang_CS5,
            py::arg("dist_mm"), R"(
            Covert the CS5 distance from the origin in mm to an angle.

            This uses the radial platescale to convert the distance in mm from the
            CS5 origin into an angle in radians.

            Args:
                dist_mm (float): The distance in mm.

            Returns:
                (float): the angle in radians.

        )")
        .def("radial_ang2dist_curved", &fba::Hardware::radial_ang2dist_curved,
            py::arg("theta_rad"), R"(
            Covert the angle from the origin to curved focal surface distance in mm.

            This uses the model S(R) arc length to convert the angle to mm in the
            curved focal surface.

            Args:
                theta_rad (float): Theta angle in radians.

            Returns:
                (float): the arc length distance in mm.

        )")
        .def("radial_dist2ang_curved", &fba::Hardware::radial_dist2ang_curved,
            py::arg("arc_mm"), R"(
            Covert the curved focal surface distance from the origin in mm to an angle.

            This uses the model S(R) arc length to convert distance in mm from the
            origin into an angle in radians.

            Args:
                arc_mm (float): The arc length distance in mm.

            Returns:
                (float): the angle in radians.

        )")
        .def("radec2xy", &fba::Hardware::radec2xy, py::arg("tilera"),
            py::arg("tiledec"), py::arg("tiletheta"), py::arg("ra"),
            py::arg("dec"), py::arg("use_CS5"), R"(
            Project an RA/DEC value into X/Y coordinates.

            For the tile pointed at (tilera, tiledec), project the (ra, dec)
            value into X/Y mm.

            Args:
                tilera (float): Tile RA
                tiledec (float): Tile DEC
                tiletheta (float): The field rotation of the tile.
                ra (float): RA to project.
                dec (float): DEC to project.
                use_CS5 (bool):  If True, use CS5 coordinates, else curved.

            Returns:
                (tuple): the (X, Y) projected location.

        )")
        .def("radec2xy_multi", [](fba::Hardware & self, double tile_ra,
                double tile_dec, double tile_theta,
                std::vector <double> const & ra,
                std::vector <double> const & dec, bool use_CS5, int threads) {
                std::vector <std::pair <double, double> > xy;
                self.radec2xy_multi(
                    tile_ra, tile_dec, tile_theta, ra, dec, xy, use_CS5, threads
                );
                return xy;
            }, py::return_value_policy::take_ownership, py::arg("tilera"),
            py::arg("tiledec"), py::arg("tiletheta"), py::arg("ra"),
            py::arg("dec"), py::arg("use_CS5"), py::arg("threads"), R"(
            Project multiple RA/DEC values into X/Y coordinates.

            For the tile pointed at (tilera, tiledec), project the (ra, dec)
            values into X/Y mm.

            Args:
                tilera (float): Tile RA
                tiledec (float): Tile DEC
                tiletheta (float): The field rotation of the tile.
                ra (array): Array of RA values to project.
                dec (float): Array of DEC values to project.
                use_CS5 (bool):  If True, use CS5 coordinates, else curved.
                threads (int): If <= 0 use maximum threads,
                    else use this number.

            Returns:
                (list): list of (X, Y) tuples with projected locations.

        )")
        .def("xy2radec", &fba::Hardware::xy2radec, py::arg("tilera"),
            py::arg("tiledec"), py::arg("tiletheta"), py::arg("x"),
            py::arg("y"), py::arg("use_CS5"), R"(
            Compute the RA/DEC value of the specified X/Y location.

            For the tile pointed at (tilera, tiledec), compute the RA/DEC
            pointing of the specified X/Y location in millimeters.

            Args:
                tilera (float): Tile RA
                tiledec (float): Tile DEC
                tiletheta (float): The field rotation of the tile.
                x (float): X position in mm.
                y (float): Y position in mm.
                use_CS5 (bool):  If True, use CS5 coordinates, else curved.

            Returns:
                (tuple): the (RA, DEC) of the focalplane location.

        )")
        .def("xy2radec_multi", [](fba::Hardware & self, double tile_ra,
                double tile_dec, double tile_theta,
                std::vector <double> const & x_mm,
                std::vector <double> const & y_mm, bool use_CS5, int threads) {
                std::vector <std::pair <double, double> > radec;
                self.xy2radec_multi(
                    tile_ra, tile_dec, tile_theta, x_mm, y_mm, radec,
                    use_CS5, threads
                );
                return radec;
            }, py::return_value_policy::take_ownership, py::arg("tilera"),
            py::arg("tiledec"), py::arg("tiletheta"), py::arg("x"),
            py::arg("y"), py::arg("use_CS5"), py::arg("threads"), R"(
            Compute the RA/DEC values of the specified X/Y locations.

            For the tile pointed at (tilera, tiledec), compute the RA/DEC
            pointings of the specified X/Y locations in millimeters.

            Args:
                tilera (float): Tile RA
                tiledec (float): Tile DEC
                tiletheta (float): The field rotation of the tile.
                x (array): Array of X positions in mm.
                y (array): Array of Y positions in mm.
                use_CS5 (bool):  If True, use CS5 coordinates, else curved.
                threads (int): If <= 0 use maximum threads,
                    else use this number.

            Returns:
                (list): list of (RA, DEC) tuples.

        )")
        .def("xy_to_thetaphi", [](
                fba::Hardware & self,
                std::pair <double, double> const & center,
                std::pair <double, double> const & xy,
                double theta_arm,
                double phi_arm,
                double theta_zero,
                double phi_zero,
                double theta_min,
                double phi_min,
                double theta_max,
                double phi_max
            ) {
                double theta;
                double phi;
                bool failed = self.xy_to_thetaphi(
                    theta, phi, center, xy, theta_arm, phi_arm, theta_zero,
                    phi_zero, theta_min, phi_min, theta_max, phi_max
                );
                if (failed) {
                    return py::make_tuple(py::none(), py::none());
                } else {
                    return py::make_tuple(theta, phi);
                }
            }, py::return_value_policy::take_ownership, py::arg("center"),
            py::arg("xy"), py::arg("theta_arm"), py::arg("phi_arm"),
            py::arg("theta_zero"), py::arg("phi_zero"),
            py::arg("theta_min"), py::arg("phi_min"),
            py::arg("theta_max"), py::arg("phi_max"), R"(
            Compute the theta / phi arm angles when moved to an x / y point.

            Note that all X/Y calculations involving positioners are performed
            in the curved focal surface (NOT CS5).

            Args:
                center (tuple): The (X, Y) tuple of the device center.
                xy (tuple): The (X, Y) tuple at which to place the fiber.
                theta_arm (float): The length of the theta arm.
                phi_arm (float): The length of the phi arm.
                theta_zero (float): The theta offset.
                phi_zero (float): The phi offset.
                theta_min (float): The theta min relative to the offset.
                phi_min (float): The phi min relative to the offset.
                theta_max (float): The theta max relative to the offset.
                phi_max (float): The phi max relative to the offset.

            Returns:
                (tuple): The theta / phi angles or (None, None) if the
                    x/y location is not reachable.

        )")
        .def("loc_position_xy", &fba::Hardware::loc_position_xy, py::arg("id"),
            py::arg("xy"), py::arg("shptheta"), py::arg("shpphi"), R"(
            Move a positioner to a given location.

            This takes the specified location and computes the shapes of
            the central body and fiber holder when the fiber is moved to
            a given (X, Y) position in the curved focal surface.  The input shapes
            are modified in place.

            Args:
                loc (int): Device location.
                xy (tuple): The (X, Y) tuple at which to place the fiber.
                shptheta (Shape):  The theta shape.
                shpphi (Shape):  The phi shape.

            Returns:
                None

        )")
        .def("loc_position_xy_multi", &fba::Hardware::loc_position_xy_multi,
            py::arg("loc"), py::arg("xy"), py::arg("threads"), R"(
            Move positioners to given locations.

            This takes the specified locations and computes the shapes of
            the central body and fiber holder when each fiber is moved to
            a given (X, Y) position in the curved focal surface.  The returned
            list of tuples contains the (central body, fiber holder) as 2 Shape
            objects for each location.

            Args:
                loc (list): List of locations.
                xy (list): List of (X, Y) tuples at which to place each fiber.
                threads (int): If <= 0 use maximum threads,
                    else use this number.

            Returns:
                (list): One tuple for each location with positioner shapes.

        )")
        .def("loc_position_thetaphi", &fba::Hardware::loc_position_thetaphi,
            py::arg("loc"), py::arg("theta"), py::arg("phi"),
            py::arg("shptheta"), py::arg("shpphi"), R"(
            Move a positioner to a given set of theta / phi angles.

            This takes the specified angles and computes the shapes of
            the central body and fiber holder when the fiber is moved to
            a the theta / phi orientation.  The input shapes are modified
            in place.

            Args:
                loc (int): Device location.
                theta (float): The theta angle.
                phi (float): The phi angle.
                shptheta (Shape):  The theta shape.
                shpphi (Shape):  The phi shape.

            Returns:
                None

        )")
        .def("check_collisions_xy", &fba::Hardware::check_collisions_xy,
            py::arg("loc"), py::arg("xy"), py::arg("threads"), R"(
            Check for collisions.

            This takes the specified locations in the curved focal surface
            and computes the shapes of the central body and fiber holder
            when each fiber is moved to a given (X, Y) position.  It then
            tests for collisions between these shapes among the locations
            specified.  The returned list of bools is True whenever the
            corresponding fiber had a collision and False otherwise.

            Args:
                loc (list): List of locations.
                xy (list): List of (X, Y) tuples at which to place each fiber.
                threads (int): If <= 0 use maximum threads,
                    else use this number.

            Returns:
                (list): A boolean value for each location.

        )")
        .def("check_collisions_thetaphi",
             &fba::Hardware::check_collisions_thetaphi, py::arg("loc"), py::arg("theta"), py::arg("phi"), py::arg("threads"), R"(
             Check for collisions.

             This takes the specified locations and computes the shapes of
             the central body and fiber holder when each fiber is moved to
             a given (theta, phi) orientation.  It then tests for collisions
             between these shapes among the locations specified.  The returned
             list of bools is True whenever the corresponding fiber had a
             collision and False otherwise.

             Args:
                 loc (array): List of locations.
                 theta (array): Theta angle for each positioner.
                 phi (array): Phi angle for each positioner.
                 threads (int): If <= 0 use maximum threads,
                     else use this number.

             Returns:
                 (list): A boolean value for each location.

         )")
        .def(py::pickle(
            [](fba::Hardware const & p) { // __getstate__
                int32_t nloc = p.locations.size();
                auto ps_radius = p.platescale_radius_mm();
                auto ps_theta = p.platescale_theta_deg();
                auto arclen = p.radial_arclen();
                std::string timestr = p.time();
                std::vector <int32_t> lid(nloc);
                std::vector <int32_t> petal(nloc);
                std::vector <int32_t> device(nloc);
                std::vector <int32_t> slitblock(nloc);
                std::vector <int32_t> blockfiber(nloc);
                std::vector <int32_t> fiber(nloc);
                std::vector <std::string> device_type(nloc);
                std::vector <double> x_mm(nloc);
                std::vector <double> y_mm(nloc);
                std::vector <int32_t> status(nloc);
                std::vector <double> theta_offset(nloc);
                std::vector <double> theta_min(nloc);
                std::vector <double> theta_max(nloc);
                std::vector <double> theta_arm(nloc);
                std::vector <double> phi_offset(nloc);
                std::vector <double> phi_min(nloc);
                std::vector <double> phi_max(nloc);
                std::vector <double> phi_arm(nloc);
                std::vector <fbg::shape> excl_theta(nloc);
                std::vector <fbg::shape> excl_phi(nloc);
                std::vector <fbg::shape> excl_gfa(nloc);
                std::vector <fbg::shape> excl_petal(nloc);
                for (int32_t i = 0; i < nloc; ++i) {
                    lid[i] = p.locations.at(i);
                    petal[i] = p.loc_petal.at(lid[i]);
                    device[i] = p.loc_device.at(lid[i]);
                    slitblock[i] = p.loc_slitblock.at(lid[i]);
                    blockfiber[i] = p.loc_blockfiber.at(lid[i]);
                    fiber[i] = p.loc_fiber.at(lid[i]);
                    device_type[i] = p.loc_device_type.at(lid[i]);
                    x_mm[i] = p.loc_pos_cs5_mm.at(lid[i]).first;
                    y_mm[i] = p.loc_pos_cs5_mm.at(lid[i]).second;
                    status[i] = p.state.at(lid[i]);
                    theta_offset[i] = p.loc_theta_offset.at(lid[i]);
                    theta_min[i] = p.loc_theta_min.at(lid[i]);
                    theta_max[i] = p.loc_theta_max.at(lid[i]);
                    theta_arm[i] = p.loc_theta_arm.at(lid[i]);
                    phi_offset[i] = p.loc_phi_offset.at(lid[i]);
                    phi_min[i] = p.loc_phi_min.at(lid[i]);
                    phi_max[i] = p.loc_phi_max.at(lid[i]);
                    phi_arm[i] = p.loc_phi_arm.at(lid[i]);
                    excl_theta[i] = p.loc_theta_excl.at(lid[i]);
                    excl_phi[i] = p.loc_phi_excl.at(lid[i]);
                    excl_gfa[i] = p.loc_gfa_excl.at(lid[i]);
                    excl_petal[i] = p.loc_petal_excl.at(lid[i]);
                }
                return py::make_tuple(
                    timestr,
                    lid, petal, device, slitblock, blockfiber, fiber,
                    device_type, x_mm, y_mm, status, theta_offset,
                    theta_min, theta_max, theta_arm, phi_offset, phi_min,
                    phi_max, phi_arm, ps_radius, ps_theta, arclen,
                    excl_theta, excl_phi,
                    excl_gfa, excl_petal);
            },
            [](py::tuple t) { // __setstate__
                return new fba::Hardware(
                    t[0].cast<std::string>(),
                    t[1].cast<std::vector<int32_t> >(),
                    t[2].cast<std::vector<int32_t> >(),
                    t[3].cast<std::vector<int32_t> >(),
                    t[4].cast<std::vector<int32_t> >(),
                    t[5].cast<std::vector<int32_t> >(),
                    t[6].cast<std::vector<int32_t> >(),
                    t[7].cast<std::vector<std::string> >(),
                    t[8].cast<std::vector<double> >(),
                    t[9].cast<std::vector<double> >(),
                    t[10].cast<std::vector<int32_t> >(),
                    t[11].cast<std::vector<double> >(),
                    t[12].cast<std::vector<double> >(),
                    t[13].cast<std::vector<double> >(),
                    t[14].cast<std::vector<double> >(),
                    t[15].cast<std::vector<double> >(),
                    t[16].cast<std::vector<double> >(),
                    t[17].cast<std::vector<double> >(),
                    t[18].cast<std::vector<double> >(),
                    t[19].cast<std::vector<double> >(),
                    t[20].cast<std::vector<double> >(),
                    t[21].cast<std::vector<double> >(),
                    t[22].cast<std::vector<fbg::shape> >(),
                    t[23].cast<std::vector<fbg::shape> >(),
                    t[24].cast<std::vector<fbg::shape> >(),
                    t[25].cast<std::vector<fbg::shape> >()
                );
            }
        ));


    py::class_ <fba::Target, fba::Target::pshr > (m, "Target", R"(
        Class representing a single target.
        )")
        .def(py::init <> ())
        .def_readonly("id", &fba::Target::id, R"(
            The target ID.
        )")
        .def_readwrite("ra", &fba::Target::ra, R"(
            The target RA.
        )")
        .def_readwrite("dec", &fba::Target::dec, R"(
            The target DEC.
        )")
        .def_readwrite("bits", &fba::Target::bits, R"(
            The target bitfield (e.g. DESI_TARGET, CMX_TARGET, etc).
        )")
        .def_readwrite("obsremain", &fba::Target::obsremain, R"(
            The remaining observations for this target.
        )")
        .def_readwrite("priority", &fba::Target::priority, R"(
            The integer priority class for this target.
        )")
        .def_readwrite("subpriority", &fba::Target::subpriority, R"(
            The float64 subpriority on the range [0,1).
        )")
        .def_readwrite("obscond", &fba::Target::obscond, R"(
            The valid observing conditions allowed for this target.
        )")
        .def_readwrite("type", &fba::Target::type, R"(
            The internal target type (science, standard, sky, safe).
        )")
        .def("is_science", &fba::Target::is_science, R"(
            Returns True if this is a science target, else False.
        )")
        .def("is_standard", &fba::Target::is_standard, R"(
            Returns True if this is a standard target, else False.
        )")
        .def("is_sky", &fba::Target::is_sky, R"(
            Returns True if this is a sky target, else False.
        )")
        .def("is_suppsky", &fba::Target::is_suppsky, R"(
            Returns True if this is a suppsky target, else False.
        )")
        .def("is_safe", &fba::Target::is_safe, R"(
            Returns True if this is a safe target, else False.
        )")
        .def("total_priority", &fba::Target::total_priority, R"(
            Return the total priority based on PRIORITY, SUBPRIORITY, and obs remaining.
        )")
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
    PYBIND11_NUMPY_DTYPE(fba::Target, id, ra, dec, bits, obsremain, priority,
                         subpriority, obscond, type);


    py::class_ <fba::Targets, fba::Targets::pshr > (m, "Targets", R"(
        Class representing a list of targets.

        The target data is stored internally as a dictionary of Target
        instances.  After construction this container is empty, and then one
        uses the "append" method to add Targets.

        )")
        .def(py::init < > ())
        .def("append", &fba::Targets::append, py::arg("tsurvey"),
            py::arg("ids"), py::arg("ras"), py::arg("decs"),
            py::arg("targetbits"), py::arg("obsremain"),
            py::arg("priority"), py::arg("subpriority"),
            py::arg("obscond"), py::arg("type"), R"(
            Append objects to the target list.

            Args:
                survey (str):  the survey type of the target data.
                ids (array):  array of int64 target IDs.
                ras (array):  array of float64 target RA coordinates.
                decs (array):  array of float64 target DEC coordinates.
                obsremain (array):  array of int32 number of remaining
                    observations.
                targetbits (array):  array of int64 bit values (DESI_TARGET,
                    CMX_TARGET, etc).
                priority (array):  array of int32 values representing the
                    target class priority for each object.
                subpriority (array):  array of float64 values in [0.0, 1.0]
                    representing the priority within the target class.
                obscond (array):  array of int32 bitfields describing the
                    valid observing conditions for each target.
                type (array):  array of uint8 bitfields holding the types of
                    of each target (science, standard, etc).
                survey (list):  list of strings of the survey types for each
                    target.

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
            }, R"(
                Returns an array of all target IDs.
        )")
        .def("get", [](fba::Targets & self, int64_t id) {
                return self.data[id];
            }, py::return_value_policy::reference_internal, py::arg("id"), R"(
            Get the specified target.

            Return a reference to the internal Target object with the
            specified ID.

            Args:
                id (int): The target ID

            Returns:
                (Target): A target object.

        )")
        .def("survey", [](fba::Targets & self) {
                return self.survey;
            }, py::return_value_policy::copy, R"(
            Get the survey for the Targets object.

            Returns:
                (str): The survey name.

        )")
        .def("__repr__",
            [](fba::Targets const & tgs) {
                std::ostringstream o;
                o << "<fiberassign.Targets with " << tgs.data.size() << " objects>";
                return o.str();
            }
        );


    py::class_ <fba::Tiles, fba::Tiles::pshr > (m, "Tiles", R"(
        Class representing a sequence of tiles.

        This is just a container of tile properties and mapping from
        tile ID to the order in the list.

        Args:
            ids (array):  array of int32 tile IDs.
            ras (array):  array of float64 tile RA coordinates.
            decs (array):  array of float64 tile DEC coordinates.
            obs (array):  array of int32 obsconditions bitfields.
            timeobs (list):  List of iso format datetime strings.
            thetaobs (array):  Array of field rotation angles.
            hourangobs (array):  Array of hour angles.

        )")
        .def(
            py::init(
                [](
                    std::vector <int32_t> const & ids,
                    std::vector <double> const & ras,
                    std::vector <double> const & decs,
                    std::vector <int32_t> const & obs,
                    py::list timeobs,
                    std::vector <double> const & thetaobs,
                    std::vector <double> const & hourangobs
                ) {
                    std::vector <std::string> tobs;
                    for (auto ts : timeobs) {
                        tobs.push_back(ts.cast <std::string> ());
                    }
                    return new fba::Tiles(
                        ids, ras, decs, obs, tobs, thetaobs, hourangobs
                    );
                }
            )
        )
        .def_readonly("id", &fba::Tiles::id, R"(
            The array of tile IDs.
        )")
        .def_readwrite("ra", &fba::Tiles::ra, R"(
            The array of tile RA values.
        )")
        .def_readwrite("dec", &fba::Tiles::dec, R"(
            The array of tile DEC values.
        )")
        .def_readwrite("obscond", &fba::Tiles::obscond, R"(
            The array of tile observing conditions values.
        )")
        .def_readwrite("obstime", &fba::Tiles::obstime, R"(
            The array of tile observation times.
        )")
        .def_readwrite("obstheta", &fba::Tiles::obstheta, R"(
            The array of tile field rotation values.
        )")
        .def_readwrite("obshourang", &fba::Tiles::obshourang, R"(
            The array of tile observation hour angles.
        )")
        .def_readonly("order", &fba::Tiles::order, R"(
            Dictionary of tile index for each tile ID.
        )")
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
            py::arg("tgs"), py::arg("min_tree_size") = 0.01)
        .def("near", [](fba::TargetTree & self, double ra_deg,
                double dec_deg, double radius_rad) {
                std::vector <int64_t> result;
                self.near(ra_deg, dec_deg, radius_rad, result);
                return result;
            }, py::return_value_policy::take_ownership, py::arg("ra"),
            py::arg("dec"), py::arg("radius"), R"(
            Get target IDs within a radius of a given point.

            Returns an array of target IDs located within the specified
            radius of the given RA/DEC.

            Args:
                ra (float): The RA of the sky location in degrees.
                dec (float): The DEC of the sky location in degrees.
                radius (float): The radius in **radians**.

            Returns:
                (array): An array of target IDs.

        )");

    py::class_ <fba::TargetsAvailable, fba::TargetsAvailable::pshr > (m,
        "TargetsAvailable", R"(
        Class representing the objects reachable by each location of each tile.

        This data structure makes it convenient and efficient to get the list
        of target IDs available to a given location for a given tile.

        Args:
            hw (Hardware):  The hardware model.
            objs (Targets):  The Targets.
            tiles (Tiles):  The tiles to consider.
            tree (TargetTree):  The HTM tree of object positions.

        )")
        .def(py::init < fba::Hardware::pshr, fba::Targets::pshr,
             fba::Tiles::pshr, fba::TargetTree::pshr > (), py::arg("hw"),
             py::arg("objs"), py::arg("tiles"), py::arg("tree")
        )
        .def("hardware", &fba::TargetsAvailable::hardware, R"(
            Return a handle to the Hardware object used.
        )")
        .def("tiles", &fba::TargetsAvailable::tiles, R"(
            Return a handle to the Tiles object used.
        )")
        .def("tile_data", &fba::TargetsAvailable::tile_data,
            py::return_value_policy::reference_internal, py::arg("tile"), R"(
            Return the targets available for a given tile.

            This returns a copy of the internal C++ available targets for a
            given tile ID.  The returned data is a dictionary with the
            location as the key and the value is an array of target IDs
            available to the location.

            Args:
                tile (int): The tile ID.

            Returns:
                (dict): Dictionary of available targets for each location.

        )");


    py::class_ <fba::LocationsAvailable, fba::LocationsAvailable::pshr > (m,
        "LocationsAvailable", R"(
        Class representing the tile/location reachable by each target.

        This data structure makes it convenient and efficient to get the list
        of tile / location pairs that can reach a given target.

        Args:
            tgsavail (TargetsAvailable):  the targets available to each
                location.

        )")
        .def(py::init < fba::TargetsAvailable::pshr > (), py::arg("tgsavail"))
        .def("target_data", &fba::LocationsAvailable::target_data,
            py::arg("target"), R"(
            Return the tile/loc pairs that can reach a target.

            This returns a copy of the internal C++ data.  The return value
            is a list of tuples containing the (tile ID, location) pairs that
            can reach the specified target ID.

            Args:
                target (int): The target ID.

            Returns:
                (list): List of (tile, loc) tuples.

        )");


    py::class_ <fba::Assignment, fba::Assignment::pshr > (m,
        "Assignment", R"(
        Class representing the current assignment of all locations.

        This data structure stores the current location assignment and
        provides methods for manipulating that assignment

        Args:
            tgs (Targets):  the targets.
            tgsavail (TargetsAvailable):  the targets available to each
                location.
            locavail (LocationsAvailable):  the locations available to each
                target.

        )")
        .def(py::init < fba::Targets::pshr, fba::TargetsAvailable::pshr,
             fba::LocationsAvailable::pshr > (), py::arg("tgs"),
             py::arg("tgsavail"), py::arg("locavail")
         )
        .def("targets", &fba::Assignment::targets, R"(
            Return a handle to the Targets object used.
        )")
        .def("hardware", &fba::Assignment::hardware, R"(
            Return a handle to the Hardware object used.
        )")
        .def("tiles", &fba::Assignment::tiles, R"(
            Return a handle to the Tiles object used.
        )")
        .def("targets_avail", &fba::Assignment::targets_avail, R"(
            Return a handle to the TargetsAvailable object used.
        )")
        .def("locations_avail", &fba::Assignment::locations_avail, R"(
            Return a handle to the LocationsAvailable object used.
        )")
        .def("tiles_assigned", &fba::Assignment::tiles_assigned, R"(
            Return an array of currently assigned tile IDs.
        )")
        .def("tile_location_target", &fba::Assignment::tile_location_target,
            py::return_value_policy::reference_internal, py::arg("tile"), R"(
            Return the assignment for a given tile.

            This returns a copy of the internal C++ target assignment for a
            given tile ID.  The returned data is a dictionary with the
            location as the key and the value is the assigned target ID.
            **Only assigned locations are stored**.  Unassigned locations do
            not exist in this dictionary.

            Args:
                tile (int): The tile ID.

            Returns:
                (dict): Dictionary of assigned target for each location.

        )")
        .def("assign_unused", &fba::Assignment::assign_unused,
             py::arg("tgtype")=TARGET_TYPE_SCIENCE,
             py::arg("max_per_petal")=-1,
             py::arg("pos_type")=std::string("POS"),
             py::arg("start_tile")=-1, py::arg("stop_tile")=-1, R"(
            Assign targets to unused locations.

            This will attempt to assign targets of the specified type to
            unused locations with devices of the specified type.

            Args:
                tgtype (int): The target type to assign, which must be one
                    of the predefined TARGET_TYPE_* module attributes.
                max_per_petal (int): Limit the assignment to this many objects
                    per petal.  Default is no limit.
                pos_type (str): Only consider this positioner device type.
                    Default is "POS".
                start_tile (int): Start assignment at this tile ID in the
                    sequence of tiles.
                stop_tile (int): Stop assignment at this tile ID (inclusive)
                    in the sequence of tiles.

            Returns:
                None

        )")
        .def("assign_force", &fba::Assignment::assign_force,
             py::arg("tgtype")=TARGET_TYPE_SCIENCE,
             py::arg("required_per_petal")=0,
             py::arg("start_tile")=-1, py::arg("stop_tile")=-1, R"(
            Force assignment of targets to unused locations.

            This function will "bump" science targets (starting with lowest
            priority) in order to place the required number of targets of the
            specified type on each petal.

            Args:
                tgtype (int): The target type to assign, which must be one
                    of the predefined TARGET_TYPE_* module attributes.
                required_per_petal (int): Bump science targets until this
                    limit is reached.
                start_tile (int): Start assignment at this tile ID in the
                    sequence of tiles.
                stop_tile (int): Stop assignment at this tile ID (inclusive)
                    in the sequence of tiles.

            Returns:
                None

        )")
        .def("redistribute_science", &fba::Assignment::redistribute_science,
             py::arg("start_tile")=-1, py::arg("stop_tile")=-1, R"(
            Redistribute science targets to future tiles.

            This function attempts to load balance the science targets per
            petal by moving science target to future available tile/loc
            placements that lie on petals with fewer total science targets.

            Args:
                start_tile (int): Start assignment at this tile ID in the
                    sequence of tiles.
                stop_tile (int): Stop assignment at this tile ID (inclusive)
                    in the sequence of tiles.

            Returns:
                None

        )");


}
