
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "fiberassign.h"


namespace py = pybind11;

using ShapeContainer = py::detail::any_container<ssize_t>;


PYBIND11_MODULE(_internal, m) {
    m.doc() = "Internal wrapper around compiled fiberassign code.";

    // As a starting point, simply define a single function to replace
    // the old fiberassign_exec main().  This can then be progressively
    // refined.

    m.def("fiberassign_exec", &fiberassign_exec);

}
