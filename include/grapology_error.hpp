// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#define np_debug_info __FILE__, __FUNCTION__, __LINE__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef GRAPOLOGYERRORHPP
#define GRAPOLOGYERRORHPP

namespace grapology {

template <typename... T>
inline void error(const T&... errs) {
  (::std::cerr << ... << errs) << "\n";
  throw std::runtime_error("Runtime error occurred");
}

template <typename... T>
inline void error(const char* FILE,
                  const char* FUNCTION,
                  const decltype(__LINE__)& LINE,
                  const T&... errs) {
  ::std::cerr << "\033[1m" << "\033[31m"
              << std::format(
                     "\nAn error was reported in the source code at 'line {2}' "
                     "within the function '{1}()' in file '{0}'\n",
                     FILE, FUNCTION, LINE)
              << "\033[0m" << "\033[0m";
  (::std::cerr << ... << errs) << "\n";
  throw std::runtime_error("Runtime error occurred");
}

}  // namespace grapology

#endif