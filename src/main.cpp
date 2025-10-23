// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: 2023-2025 Zhengrong Guo

#include <chrono>
#include "grapology.h"  // include this to use the library

int main(int argc, char* argv[]) {
  const auto start_time{std::chrono::steady_clock::now()};
  size_t line_count{0}, success_count{0};
  std::cerr << "Graphene-topology: A topological vector library for "
               "generating graphene-derived nanostructures\n";

  if (argc >= 2) {
    bool is_mol2 = false;
    bool is_xyz = false;
    bool is_pdb = false;
    bool is_lammps = false;
    bool is_lammpstrj = false;
    bool is_poscar = false;
    bool is_povray = false;
    bool is_none = false;
    bool whether_pov2 = false;
    bool whether_optimize = false;
    bool is_optimize_in_formation = true;

    for (int i = 2; i < argc; i++) {
      if (std::string_view(argv[i]) == "-mol2")
        is_mol2 = true;
      else if (std::string_view(argv[i]) == "-xyz")
        is_xyz = true;
      else if (std::string_view(argv[i]) == "-lammps")
        is_lammps = true;
      else if (std::string_view(argv[i]) == "-lammpstrj")
        is_lammpstrj = true;
      else if (std::string_view(argv[i]) == "-pdb")
        is_pdb = true;
      else if (std::string_view(argv[i]) == "-povray")
        is_povray = true;
      else if (std::string_view(argv[i]) == "-poscar")
        is_poscar = true;
      else if (std::string_view(argv[i]) == "-optimize")
        whether_optimize = true;
      else if (std::string_view(argv[i]) == "-pov2")
        whether_pov2 = true;  // for test
      else if (std::string_view(argv[i]) == "-none") {
        is_none = true;
      } else {
        std::cerr << std::format("Unknown option argument: {}\n", argv[i]);
        std::cerr << "!  Supported options: -mol2, -xyz, -pdb, -lammps, "
                     "-lammpstrj -poscar, -povray -optimize\n";
        return 1;
      }
    }

    if (is_mol2 == false and is_xyz == false and is_pdb == false and
        is_lammps == false and is_lammpstrj == false and is_poscar == false and
        is_povray == false and is_none == false)
      is_mol2 = true;
    if (is_none) {
      is_mol2 = false;
      is_xyz = false;
      is_pdb = false;
      is_lammps = false;
      is_lammpstrj = false;
      is_poscar = false;
      is_povray = false;
      whether_optimize = false;
      is_optimize_in_formation = false;
    }

    std::fstream in;
    in.open(argv[1], std::ios::in);
    if (in.is_open() == false) {
      std::cerr << std::format("Failed to open file: {}\n", argv[1]);
      return 0;
    }

    // create a carbon cap object
    grapology::carbon_cap cap;
    cap.whether_optimizing = is_optimize_in_formation;
    // create a locate object and pass a reference of atoms to it
    grapology::locate<false> loc(cap.get_atoms());
    std::string line;
    std::vector<int> paras;
    paras.reserve(13);
    std::vector<std::string> results;
    results.reserve(50);
    while (std::getline(in, line) and
           in) {  // start loop over each line in the file
      line_count++;

      std::istringstream iss(line);
      results.clear();
      for (std::string s; iss >> s;)
        results.push_back(s);
      if (results.size() < 13) {
        std::cerr << std::format("Formatting error in line: {}\n", line);
        continue;
      }

      try {
        paras.clear();
        for (size_t i = 0; i < 13; i++) {
          paras.push_back(std::stoi(results[i]));
        }
      } catch (...) {
        std::cerr << std::format("Formatting error in line {}: {}\n",
                                 line_count, line);
        std::cerr << "Please check the format of the line.\n";
        std::cerr << "!  Each line should contain 13 integers, where the first "
                     "12 integers\n"
                     "!  represent the chiral indices of carbon caps, and the "
                     "13th integer\n"
                     "!  indicates the number of additional hexagonal carbon "
                     "rings to be added\n"
                     "!  to the edge of the cap.\n";
        continue;
      }
      // clear the atoms in the cap object before constructing a new one
      cap.clear();

      try {
        // construct the carbon cap based on the specifications in the line
        cap.construct(paras);
        if (whether_optimize)
          loc.dynamic_relaxing(10000);  // optimize the structure if required
        loc.set_toward_z();      // set the cap to be oriented along the z-axis
        loc.set_edge_to_zero();  // set the edge atoms to be at the origin
      } catch (std::exception& e) {
        std::cerr << std::format(
            "Failed to generate the carbon cap based on the specifications in "
            "line {}: '{}'\n",
            line_count, line);
        continue;
      }
      success_count++;

      // save the generated carbon nanotube based on the specified format
      // given by options controling args or default setting
      if (is_mol2)
        cap.save_atoms_by_mol2(std::format("cnt.{}.mol2", success_count));
      if (is_xyz)
        cap.save_atoms_by_xyz(std::format("cnt.{}.xyz", success_count));
      if (is_pdb)
        cap.save_atoms_by_pdb(std::format("cnt.{}.pdb", success_count));
      if (is_lammps)
        cap.save_atoms_by_lammps(std::format("cnt.{}.lammps", success_count));
      if (is_lammpstrj)
        cap.save_atoms_by_lammpstrj(
            std::format("cnt.{}.lammpstrj", success_count));
      if (is_poscar)
        cap.save_atoms_by_poscar(std::format("cnt.{}.poscar", success_count));
      if (is_povray)
        cap.save_atoms_by_povray(std::format("cnt.{}.pov", success_count));
      if (whether_pov2)
        cap.save_atoms_and_edge_by_povray(
            std::format("cnt.{}.pov", success_count));
    }  // end loop over each line in the file
    in.close();
  } else {
    std::cerr
        << "No input script found! Graphene-topology reads the chiral indices "
           "from a script file. For example: ./grapology caps.txt -xyz"
        << std::endl;
    return 0;
  }

  const auto complete_time{std::chrono::steady_clock::now()};
  const std::chrono::duration<float> elapsed_seconds{complete_time -
                                                     start_time};
  std::cerr << std::format(
                   "\nGenerated {} structure{}, with a time cost of {}.",
                   success_count, success_count <= 1 ? "" : "s",
                   elapsed_seconds)
            << std::endl;
  return 0;
}