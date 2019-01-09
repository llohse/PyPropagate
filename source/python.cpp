#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "finite_difference.h"

namespace py = pybind11;

PYBIND11_MODULE(_pypropagate, m){
  
  py::class_<finite_difference_AF>(m, "finite_difference_AF")
  .def_readwrite("ra",&finite_difference_AF::ra)
  .def_readwrite("rf",&finite_difference_AF::rf)
  .def_readwrite("u",&finite_difference_AF::u)
  .def(py::init<>())
  .def("step", &finite_difference_AF::step)
  .def("update", &finite_difference_AF::update)
  .def("resize", &finite_difference_AF::resize)
  ;
  
  py::class_<finite_difference_ACF>(m, "finite_difference_ACF")
  .def_readwrite("ra",&finite_difference_ACF::ra)
  .def_readwrite("rc",&finite_difference_ACF::rc)
  .def_readwrite("rf",&finite_difference_ACF::rf)
  .def_readwrite("u",&finite_difference_ACF::u)
  .def_readwrite("thread_count",&finite_difference_ACF::thread_count)
  .def(py::init<>())
  .def("step_1", &finite_difference_ACF::step_1)
  .def("step_2", &finite_difference_ACF::step_2)
  .def("update", &finite_difference_ACF::update)
  .def("resize", &finite_difference_ACF::resize)
  ;
  
  py::class_<finite_difference_A0F>(m, "finite_difference_A0F")
  .def_readwrite("ra",&finite_difference_A0F::ra)
  .def_readwrite("rf",&finite_difference_A0F::rf)
  .def_readwrite("u",&finite_difference_A0F::u)
  .def_readwrite("thread_count",&finite_difference_A0F::thread_count)
  .def(py::init<>())
  .def("step", &finite_difference_A0F::step)
  .def("update", &finite_difference_A0F::update)
  .def("resize", &finite_difference_A0F::resize)
  ;

  py::class_<finite_difference_ABC>(m, "finite_difference_ABC")
  .def_readwrite("ra",&finite_difference_ABC::ra)
  .def_readwrite("rb",&finite_difference_ABC::rb)
  .def_readwrite("rc",&finite_difference_ABC::rc)
  .def_readwrite("rz",&finite_difference_ABC::rz)
  .def_readwrite("thread_count",&finite_difference_ABC::thread_count)
  .def_readwrite("u",&finite_difference_ABC::u)
  .def(py::init<>())
  .def("step", &finite_difference_ABC::step)
  .def("update", &finite_difference_ABC::update)
  .def("resize", &finite_difference_ABC::resize)
  ;
}
