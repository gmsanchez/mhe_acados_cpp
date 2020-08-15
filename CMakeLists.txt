cmake_minimum_required(VERSION 2.8.3)
project(mhe_acados_cpp)

find_package(Eigen3 REQUIRED)
include_directories(${Eigen3_INCLUDE_DIRS})

set (CMAKE_C_FLAGS "-fPIC -DACADOS_WITH_QPOASES")
# set (CMAKE_CXX_FLAGS "-fPIC -DACADOS_WITH_QPOASES")

include_directories(
# include
  "/home/guiss/fun/acados/include/blasfeo/include/"
  "/home/guiss/fun/acados/include/hpipm/include/"
  "/home/guiss/fun/acados/include"
  "/home/guiss/fun/acados/include/qpOASES_e/"
  "/home/guiss/fun/acados/"
  "/usr/include/eigen3/"
  "solver/"
  "solver/mhe/c_generated_code/"
  ${catkin_INCLUDE_DIRS}
)

link_directories(
  "/home/guiss/fun/acados/lib"
)

## Declare a C++ library
add_library(${PROJECT_NAME}_lib
#   src/${PROJECT_NAME}/agv_nonlinear_mhe_acados.cpp
  "solver/mhe/c_generated_code/mhe_pendulum_ode_model/mhe_pendulum_ode_expl_ode_fun.c"
  "solver/mhe/c_generated_code/mhe_pendulum_ode_model/mhe_pendulum_ode_expl_vde_forw.c"
  "solver/mhe/c_generated_code/mhe_pendulum_ode_cost/mhe_pendulum_ode_cost_y_fun.c"
  "solver/mhe/c_generated_code/mhe_pendulum_ode_cost/mhe_pendulum_ode_cost_y_fun_jac_ut_xt.c"
  "solver/mhe/c_generated_code/mhe_pendulum_ode_cost/mhe_pendulum_ode_cost_y_hess.c"
  "solver/mhe/c_generated_code/acados_solver_mhe_pendulum_ode.c"
  "solver/mhe/c_generated_code/acados_sim_solver_mhe_pendulum_ode.c"
)

add_executable(${PROJECT_NAME} src/mhe_pendulum_ode_test.cpp)

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)

target_link_libraries(${PROJECT_NAME} acados hpipm blasfeo qpOASES_e m)


target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
)
