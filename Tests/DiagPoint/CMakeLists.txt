cmake_minimum_required(VERSION 3.18)
project(PeleCDiagPointTests LANGUAGES CXX)
enable_testing()
# Keep a convenient PeleC entry point for the shared PelePhysics tests.
include(${CMAKE_CURRENT_LIST_DIR}/../../Submodules/PelePhysics/Testing/Diagnostics/DiagPoint/CMakeLists.txt)
