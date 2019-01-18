#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#include "finite_difference.h"

namespace py = pybind11;

using namespace finite_differences;

PYBIND11_MODULE(_pypropagate, m){
  
  py::class_<finite_difference_AF>(m, "finite_difference_AF")
  .def(py::init<size_t>())
  .def("step", &finite_difference_AF::step)
  .def("update", &finite_difference_AF::update)
  .def("get_field", &finite_difference_AF::get_field, py::return_value_policy::reference_internal)
  .def("set_field", &finite_difference_AF::set_field)
  .def("set_ra", &finite_difference_AF::set_ra)
  .def("set_rf", &finite_difference_AF::set_rf)
  ;
  
  py::class_<finite_difference_ACF>(m, "finite_difference_ACF")
  .def_readwrite("thread_count",&finite_difference_ACF::thread_count)
  .def(py::init<size_t,size_t>())
  .def("step_1", &finite_difference_ACF::step_1)
  .def("step_2", &finite_difference_ACF::step_2)
  .def("update", &finite_difference_ACF::update)
  .def("get_field", &finite_difference_ACF::get_field, py::return_value_policy::reference_internal)
  .def("set_field", &finite_difference_ACF::set_field)
  .def("set_ra", &finite_difference_ACF::set_ra)
  .def("set_rc", &finite_difference_ACF::set_rc)
  .def("set_rf", &finite_difference_ACF::set_rf)

  ;
  
  py::class_<finite_difference_A0F>(m, "finite_difference_A0F")
  .def_readwrite("thread_count",&finite_difference_A0F::thread_count)
  .def(py::init<size_t, size_t>())
  .def("step", &finite_difference_A0F::step)
  .def("update", &finite_difference_A0F::update)
  .def("get_field", &finite_difference_A0F::get_field, py::return_value_policy::reference_internal)
  .def("set_field", &finite_difference_A0F::set_field)
  .def("set_ra", &finite_difference_A0F::set_ra)
  .def("set_rf", &finite_difference_A0F::set_rf)
  ;

  py::class_<finite_difference_ABC>(m, "finite_difference_ABC")
  .def_readwrite("thread_count",&finite_difference_ABC::thread_count)
  .def(py::init<size_t,size_t>())
  .def("step", &finite_difference_ABC::step)
  .def("update", &finite_difference_ABC::update)
  .def("get_field", &finite_difference_ABC::get_field, py::return_value_policy::reference_internal)
  .def("set_field", &finite_difference_ABC::set_field)
  .def("set_ra", &finite_difference_ABC::set_ra)
  .def("set_rb", &finite_difference_ABC::set_rb)
  .def("set_rc", &finite_difference_ABC::set_rc)
  .def("set_rz", &finite_difference_ABC::set_rz)
  ;
}
