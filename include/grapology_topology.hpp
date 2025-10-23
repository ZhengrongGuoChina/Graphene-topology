// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGY_TOPOLOGY_HPP
#define GRAPOLOGY_TOPOLOGY_HPP

namespace grapology {

/**
 * @brief Template class that implements topological representation and
 * manipulation of graphene-derived nanostructures
 *
 * This class is the core of the Grapology library, providing methods to create,
 * manipulate, and analyze topological vectors, allowing static polymorphism and
 * derived class access through template specialization.
 */
using atom_type = vertex<topo_vector>;
using atom_container = std::unordered_set<atom_type*>;
using vector_container = std::unordered_set<topo_vector*>;
using topo_path_type = circular_topo_path<topo_vector*>;

template <typename T>
struct topology {
  using derived_type = T;
  /**
   * @brief Generic interface method for constructing topological structures
   * @param paras Integer vector of construction parameters, specifics defined
   * by derived classes
   *
   * This is an implementation based on the Curiously Recurring Template
   * Pattern, which provides a unified interface in the base class while
   * delegating concrete structure building logic to derived classes.
   */
  void construct(const std::vector<int>& paras) {
    static_assert(std::is_base_of_v<topology, derived_type>,
                  "derived_type must be derived from topology<derived_type>");
    auto derived = static_cast<derived_type*>(this);
    derived->construct_derived_structure(paras);
  }

  /**
   * @brief Exports the current structure to a MOL2 file format
   * @param filename Name of the file to save the structure to
   *
   * This method saves all atoms in the structure to a file in MOL2 format. The
   * MOL2 format is a common file format for molecular structures, including
   * atomic coordinates, atom types, and bond information.
   *
   * @note The method includes a comment line "#Graphene-Derived Nanostructures"
   * in the second line of the file to identify the structure as a
   * graphene-derived nanostructure.
   */
  void save_atoms_by_mol2(const std::string& filename) const;

  /**
   * @brief Exports the current structure to a XYZ file format
   * @param filename Name of the file to save the structure to
   *
   * This method saves all atoms in the structure to a file in XYZ format. The
   * XYZ format includes the total number of atoms on the first line, a comment
   * line on the second line, and then one line per atom containing the atom
   * label and its x, y, z coordinates. The comment line contains the text
   * "#Graphene-Derived Nanostructures".
   */
  void save_atoms_by_xyz(const std::string& filename) const;

  /**
   * @brief Exports the current structure to a PDB file format
   * @param filename Name of the file to save the structure to
   *
   * This method saves all atoms in the structure to a file in PDB format. The
   * PDB format is a common file format for molecular structures, including
   * atomic coordinates, atom types, and bond information.
   *
   * @note The method includes a comment line "#Graphene-Derived Nanostructures"
   * in the second line of the file to identify the structure as a
   * graphene-derived nanostructure.
   */
  void save_atoms_by_pdb(const std::string& filename) const;

  /**
   * @brief Exports the current structure to a LAMMPS data file format or
   * lammpstrj file format
   * @param filename Name of the file to save the structure to
   *
   * This method saves all atoms in the structure to a file in LAMMPS data
   * format. The LAMMPS data format is a common file format for molecular
   * simulations, including atomic coordinates and atom types.
   *
   * @note The method includes a comment line "#Graphene-Derived Nanostructures"
   * in the second line of the file to identify the structure as a
   * graphene-derived nanostructure.
   */
  void save_atoms_by_lammps(const std::string& filename) const;
  void save_atoms_by_lammpstrj(const std::string& filename) const;

  /**
   * @brief Exports the current structure to a POV-Ray file format
   * @param filename Name of the file to save the structure to
   *
   * This method saves all atoms in the structure to a file in POV-Ray
   * (Persistence of Vision Raytracer) format. The POV-Ray format is used for
   * creating high-quality 3D renderings of molecular structures. The file is
   * quite primitive, and the user is advised to modify it according to their
   * needs. The exported file includes:
   * - Camera setup with automatic positioning based on the structure dimensions
   * - Lighting configuration
   * - Material definitions for carbon atoms and bonds
   * - Geometric primitives representing atoms (spheres) and bonds (cylinders)
   */
  void save_atoms_by_povray(const std::string& filename) const;
  // debug method
  void save_atoms_and_edge_by_povray(const std::string& filename) const;
  /**
   * @brief Saves the atomic structure to a file in POSCAR format
   *
   * This method exports the atomic coordinates and structure information
   * to a file using the POSCAR format, which is commonly used in materials
   * science software (e.g., VASP). The function is marked as const to indicate
   * it does not modify the internal state of the object.
   *
   * @param filename Name of the file to save the structure to
   */
  void save_atoms_by_poscar(const std::string& filename) const;

  /**
   * @brief Gets a constant reference to the container holding all atoms
   *
   * This method provides read-only access to the internal atom container,
   * ensuring that callers cannot modify the container contents. The function is
   * marked as constexpr and inline for compile-time evaluation and efficient
   * inlining.
   *
   * @return const atom_container& Constant reference to the atom container
   */
  constexpr inline auto get_atoms() const -> const atom_container& {
    return atoms;
  }

  /**
   * @brief Gets a constant reference to the container holding all vectors
   *
   * This method provides read-only access to the internal vector container,
   * ensuring that callers cannot modify the container contents. The function is
   * marked as constexpr and inline for compile-time evaluation and efficient
   * inlining.
   *
   * @return const vector_container& Constant reference to the vector container
   */
  constexpr inline auto get_vectors() const -> const vector_container& {
    return vectors;
  }

  /**
   * @brief Cleans up all dynamically allocated resources and empties containers
   *
   * This method performs a complete resource cleanup operation including:
   * 1. Deallocating memory for all atom objects
   * 2. Deallocating memory for all vector objects
   * 3. Clearing the contents of both atom and vector containers
   *
   * This function ensures no memory leaks occur by explicitly deleting each
   * dynamically allocated object before clearing the corresponding containers.
   */
  void clear() {
    for (auto i : atoms) {
      i->~atom_type();
      atom_pool.free(i);
    }
    for (auto v : vectors) {
      v->~topo_vector();
      vector_pool.free(v);
    }
    atoms.clear();
    vectors.clear();
  }

 protected:
  auto create_atom() -> atom_type*;
  void delete_atom(atom_type*);
  auto create_vector(topo_state) -> topo_vector*;
  void delete_vector(topo_vector*);
  auto create_heptagon_path() -> topo_path_type;
  auto create_pentagon_path() -> topo_path_type;
  auto create_forward_path_to(const int, const int) -> topo_path_type;
  auto create_backward_path_of(topo_path_type&) -> topo_path_type;
  auto remove_successive_reverse_vector(topo_vector*&) -> bool;
  auto find_most_active_site(topo_vector*) const -> topo_path_type;
  auto find_active_site_contains(topo_vector*) const -> topo_path_type;
  auto add_hexagon_onto_active_site(topo_path_type&) -> topo_vector*;
  auto dismantle_hexagons_from_edge_path(topo_vector*) -> topo_vector*;
  auto get_number_of_vectors_between(topo_vector*, topo_vector*) const
      -> size_t;
  auto get_number_on_the_edgepath(topo_vector*) const -> size_t;
  auto get_sum_from_path(topo_vector*) const -> std::array<int, 3>;
  topology()
      : atom_pool(sizeof(atom_type), 128),
        vector_pool(sizeof(topo_vector), 128) {}
  virtual ~topology() {
    for (auto i : atoms) {
      i->~atom_type();
      atom_pool.free(i);
    }
    for (auto v : vectors) {
      v->~topo_vector();
      vector_pool.free(v);
    }
  }

 private:
  atom_container atoms;
  vector_container vectors;
  grapology_pool atom_pool;
  grapology_pool vector_pool;
};  // class topology

template <typename T>
void topology<T>::save_atoms_by_mol2(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
  }

  std::set<std::array<size_t, 2>> bonds;
  for (auto atom : atoms) {
    for (auto nei : atom->conn_vertex) {
      if (atom->label < nei->label)
        bonds.insert({atom->label, nei->label});
      if (atoms.contains(nei) == false)
        error(np_debug_info, "The neighbor atom is not found!");
    }
  }

  out << "@<TRIPOS>MOLECULE\n";
  out << "generated by Graphene-topology library\n";
  out << std::format("{} {} 1 0 0\n", atoms.size(), bonds.size());
  out << "SMALL\nUSER_CHARGES\n";

  out << "\n@<TRIPOS>ATOM\n";
  for (auto atom : atoms) {
    out << std::format(
        "{:6d} C  {:10.4f} {:10.4f} {:10.4f} C.2  1 GDN {:6.4f}\n", atom->label,
        atom->x[0], atom->x[1], atom->x[2], 0.0000);
  }

  out << "\n@<TRIPOS>BOND\n";
  count = 0;
  for (auto bond : bonds) {
    out << std::format("{:6d} {:6d} {:6d} 1\n", ++count, bond[0], bond[1]);
  }

  out << "\n@<TRIPOS>SUBSTRUCTURE\n";
  out << "1 GDN\t 1 RESIDUE  1";
  out.close();
}

template <typename T>
void topology<T>::save_atoms_by_xyz(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
  }

  out << std::format("{}\n", atoms.size());
  out << "#generated by Graphene-topology library\n";
  for (auto atom : atoms) {
    out << std::format("{}  {:10.4f} {:10.4f} {:10.4f}\n", "C", atom->x[0],
                       atom->x[1], atom->x[2]);
  }
  out.close();
}

template <typename T>
void topology<T>::save_atoms_by_pdb(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
  }

  std::set<std::array<size_t, 2>> bonds;
  for (auto atom : atoms) {
    for (auto nei : atom->conn_vertex) {
      if (atom->label < nei->label)
        bonds.insert({atom->label, nei->label});
      if (atoms.contains(nei) == false)
        error(np_debug_info, "Neighbor atom not found!");
    }
  }

  for (auto atom : atoms) {
    out << std::format(
        "HETATM {:4d}  C   GDN     1    {:8.3f}{:8.3f}{:8.3f}  0.00  0.00           C  \n",
        atom->label, atom->x[0], atom->x[1], atom->x[2]);
  }

  for (auto atom : atoms) {
    auto&& neigborings = atom->get_bonds();
    if (neigborings.size() != 0) {
      out << std::format("CONECT{:5d}", atom->label);
      for (auto nei : neigborings)
        out << std::format("{:5d}", nei->label);
      out << "\n";
    }
  }

  out << "ENDMDL\n";
  out << "END\n";

  out.close();
}

template <typename T>
void topology<T>::save_atoms_by_lammps(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
  }

  out << "# LAMMPS data file (atomic) generated by Graphene-topology "
         "library\n\n";
  out << std::format("{} atoms\n", atoms.size());
  out << "1 atom types\n\n";

  DataTypeFloat min_x = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat min_y = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat min_z = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat max_x = std::numeric_limits<DataTypeFloat>::lowest();
  DataTypeFloat max_y = std::numeric_limits<DataTypeFloat>::lowest();
  DataTypeFloat max_z = std::numeric_limits<DataTypeFloat>::lowest();

  for (auto atom : atoms) {
    min_x = std::min(min_x, atom->x[0]);
    min_y = std::min(min_y, atom->x[1]);
    min_z = std::min(min_z, atom->x[2]);
    max_x = std::max(max_x, atom->x[0]);
    max_y = std::max(max_y, atom->x[1]);
    max_z = std::max(max_z, atom->x[2]);
  }

  DataTypeFloat buffer = 1.7;
  out << std::format("{:10.5f} {:10.5f} xlo xhi\n", min_x - buffer,
                     max_x + buffer);
  out << std::format("{:10.5f} {:10.5f} ylo yhi\n", min_y - buffer,
                     max_y + buffer);
  out << std::format("{:10.5f} {:10.5f} zlo zhi\n", min_z - buffer,
                     max_z + buffer);

  out << "\nMasses\n\n";
  out << "1 12.01070 # Carbon\n";

  out << "\nAtoms # atomic\n\n";
  for (auto atom : atoms) {
    out << std::format("{:5d} 1 {:15.8f} {:15.8f} {:15.8f}\n", atom->label,
                       atom->x[0], atom->x[1], atom->x[2]);
  }

  out.close();
}

template <typename T>
void topology<T>::save_atoms_by_lammpstrj(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
  }

  out << "ITEM: TIMESTEP\n";
  out << "0\n";

  out << "ITEM: NUMBER OF ATOMS\n";
  out << std::format("{}\n", atoms.size());

  DataTypeFloat min_x = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat min_y = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat min_z = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat max_x = std::numeric_limits<DataTypeFloat>::lowest();
  DataTypeFloat max_y = std::numeric_limits<DataTypeFloat>::lowest();
  DataTypeFloat max_z = std::numeric_limits<DataTypeFloat>::lowest();

  for (auto atom : atoms) {
    min_x = std::min(min_x, atom->x[0]);
    min_y = std::min(min_y, atom->x[1]);
    min_z = std::min(min_z, atom->x[2]);
    max_x = std::max(max_x, atom->x[0]);
    max_y = std::max(max_y, atom->x[1]);
    max_z = std::max(max_z, atom->x[2]);
  }

  DataTypeFloat buffer = 1.7;
  out << "ITEM: BOX BOUNDS pp pp pp\n";
  out << std::format("{:10.5f} {:10.5f}\n", min_x - buffer, max_x + buffer);
  out << std::format("{:10.5f} {:10.5f}\n", min_y - buffer, max_y + buffer);
  out << std::format("{:10.5f} {:10.5f}\n", min_z - buffer, max_z + buffer);

  // 输出原子数据
  out << "ITEM: ATOMS id type x y z\n";
  for (auto atom : atoms) {
    out << std::format("{:5d} 1 {:15.8f} {:15.8f} {:15.8f}\n", atom->label,
                       atom->x[0], atom->x[1], atom->x[2]);
  }

  out.close();
}

template <typename T>
void topology<T>::save_atoms_by_povray(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  DataTypeFloat max_r = std::numeric_limits<DataTypeFloat>::lowest();

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
    max_r = std::max(max_r, atom->x[0]);
    max_r = std::max(max_r, atom->x[1]);
  }

  std::set<std::array<atom_type*, 2>> bonds;
  for (auto atom : atoms) {
    for (auto nei : atom->conn_vertex) {
      if (atom->label < nei->label)
        bonds.insert({atom, nei});
      if (atoms.contains(nei) == false)
        error(np_debug_info, "Neighbor atom not found!");
    }
  }

  out << "// Graphene-Derived Nanostructures - Generated by Graphene-topology "
         "library"
         "library\n";
  out << "#include \"colors.inc\"\n";
  out << "#include \"finish.inc\"\n";
  out << "#include \"textures.inc\"\n\n";

  out << "camera {\n";
  out << "  angle 60\n";
  out << std::format("  location <0, 0, -{}>\n",
                     std::ceil((max_r + 0.3) * 1.7321));
  out << "  look_at <0, 0, 0>\n";
  out << "  right x*image_width/image_height\n";
  out << "}\n\n";

  out << "light_source {\n";
  out << "  <0, 5, -20>\n";
  out << "  color White\n";
  out << "  shadowless\n";
  out << "}\n\n";

  out << "background {color White}\n\n";

  out << "#declare CarbonAtomMaterial = material {\n";
  out << "  texture {\n";
  out << "    pigment { color <0.53,0.81,0.92> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare CarbonBondMaterial = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0.32, 0.32, 0.32> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare CarbonAtomSphereRadius = 0.3;\n";
  out << "#declare CarbonBondCylinderRadius = 0.1;\n\n";

  out << "// Chemical bonds\n";
  for (auto bond : bonds) {
    auto atom1 = bond[0];
    auto atom2 = bond[1];

    out << "cylinder {\n";
    out << std::format(
        "  <{}, {}, {}>, <{}, {}, {}>, CarbonBondCylinderRadius\n", atom1->x[0],
        atom1->x[1], atom1->x[2], atom2->x[0], atom2->x[1], atom2->x[2]);
    out << "  material { CarbonBondMaterial }\n";
    out << "}\n\n";
  }

  out << "// Carbon atoms\n";
  for (auto atom : atoms) {
    out << "sphere {\n";
    out << std::format("  <{}, {}, {}>, CarbonAtomSphereRadius\n", atom->x[0],
                       atom->x[1], atom->x[2]);
    out << "  material { CarbonAtomMaterial }\n";
    out << "}\n\n";
  }

  out.close();
}

/**
 * @brief Exports the current structure to a POSCAR file format
 * @param filename Name of the file to save the structure to
 *
 * This method saves all atoms in the structure to a file in POSCAR format,
 * which is commonly used by the VASP (Vienna Ab-initio Simulation Package)
 * for electronic structure calculations. The file includes:
 * - A comment line identifying the structure
 * - Scaling factor (set to 1.0)
 * - Lattice vectors automatically calculated from atom coordinates
 * - Element type (C for carbon)
 * - Number of carbon atoms
 * - Cartesian coordinates of all atoms
 */
template <typename T>
void topology<T>::save_atoms_by_poscar(const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  int count = 0;
  for (auto atom : atoms) {
    atom->label = ++count;
  }

  DataTypeFloat min_x = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat min_y = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat min_z = std::numeric_limits<DataTypeFloat>::max();
  DataTypeFloat max_x = std::numeric_limits<DataTypeFloat>::lowest();
  DataTypeFloat max_y = std::numeric_limits<DataTypeFloat>::lowest();
  DataTypeFloat max_z = std::numeric_limits<DataTypeFloat>::lowest();

  for (auto atom : atoms) {
    min_x = std::min(min_x, atom->x[0]);
    min_y = std::min(min_y, atom->x[1]);
    min_z = std::min(min_z, atom->x[2]);
    max_x = std::max(max_x, atom->x[0]);
    max_y = std::max(max_y, atom->x[1]);
    max_z = std::max(max_z, atom->x[2]);
  }

  DataTypeFloat buffer = 1.7;
  DataTypeFloat box_size_x = max_x - min_x + 2 * buffer;
  DataTypeFloat box_size_y = max_y - min_y + 2 * buffer;
  DataTypeFloat box_size_z = max_z - min_z + 2 * buffer;

  out << "POSCAR file written by Graphene-topology library\n";
  out << "1.0\n";

  out << std::format("{:15.8f} {:15.8f} {:15.8f}\n", box_size_x, 0.0, 0.0);
  out << std::format("{:15.8f} {:15.8f} {:15.8f}\n", 0.0, box_size_y, 0.0);
  out << std::format("{:15.8f} {:15.8f} {:15.8f}\n", 0.0, 0.0, box_size_z);

  out << "C\n";
  out << atoms.size() << "\n";

  out << "Cartesian\n";

  for (auto atom : atoms) {
    out << std::format("{:15.8f} {:15.8f} {:15.8f}\n",
                       atom->x[0] - min_x + buffer, atom->x[1] - min_y + buffer,
                       atom->x[2] - min_z + buffer);
  }

  out.close();
}

template <typename T>
void topology<T>::save_atoms_and_edge_by_povray(
    const std::string& filename) const {
  std::ofstream out;
  out.open(filename, std::ios::out);

  DataTypeFloat max_r = std::numeric_limits<DataTypeFloat>::lowest();
  int count = 0;
  for (auto atom : get_atoms()) {
    atom->label = ++count;
    max_r = std::max(max_r, static_cast<DataTypeFloat>(
                                std::abs(atom->x[0]) * 1.7321 - atom->x[2]));
    max_r = std::max(max_r, static_cast<DataTypeFloat>(
                                std::abs(atom->x[1]) * 1.7321 - atom->x[2]));
  }

  std::set<std::array<vertex<topo_vector>*, 2>> bonds;
  for (auto atom : get_atoms()) {
    atom->whether_located = true;
    for (auto nei : atom->conn_vertex) {
      if (atom->label < nei->label)
        bonds.insert({atom, nei});
      if (get_atoms().contains(nei) == false)
        error(np_debug_info, "Neighbor atom not found!");
    }
  }

  for (auto current : get_vectors()) {
    current->front_atom->whether_located = false;
    current->back_atom->whether_located = false;
  }

  out << "// Graphene-Derived Nanostructures - Generated by Grapology "
         "library\n";
  out << "#include \"colors.inc\"\n";
  out << "#include \"finish.inc\"\n";
  out << "#include \"textures.inc\"\n\n";

  out << "camera {\n";
  out << "  angle 60\n";
  out << std::format("  location <0, 0, -{}>\n", max_r + 1.0);
  out << "  look_at <0, 0, 0>\n";
  out << "  right x*image_width/image_height\n";
  out << "}\n\n";

  out << "light_source {\n";
  out << "  <5, 5, -20>\n";
  out << "  color White\n";
  out << "  shadowless\n";
  out << "}\n\n";

  out << "background {color White}\n\n";

  out << "#declare CarbonAtomMaterial = material {\n";
  out << "  texture {\n";
  out << "    pigment { color <0.53,0.81,0.92> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare CarbonBondMaterial = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0.32, 0.32, 0.32> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare ArrowStateI = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0.9, 0.32, 0.32> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare ArrowStateJ = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0.32, 0.9, 0.32> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare ArrowStateK = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0.32, 0.32, 0.9> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.5\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare ArrowStateNI = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0.9, 0., 0.> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.6\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare ArrowStateNJ = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0., 0.9, 0.> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.6\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  out << "#declare ArrowStateNK = material {\n";
  out << "  texture {\n";
  out << "    pigment { color rgb <0., 0., 0.9> }\n";
  out << "    finish {\n";
  out << "      phong 1\n";
  out << "      phong_size 50\n";
  out << "      diffuse 0.6\n";
  out << "    }\n";
  out << "  }\n";
  out << "}\n\n";

  const float arrow_radius = 0.05f;
  const float cylinder_radius = 0.05f;
  const float sphere_radius = 0.1f;

  out << std::format("#declare CarbonAtomSphereRadius = {};\n", sphere_radius);
  out << std::format("#declare CarbonBondCylinderRadius = {};\n",
                     cylinder_radius);
  out << std::format("#declare CarbonArrowRadius = {};\n\n", arrow_radius);

  out << "// Chemical bonds\n";
  for (auto bond : bonds) {
    auto atom1 = bond[0];
    auto atom2 = bond[1];
    if (atom1->whether_located == false and atom2->whether_located == false)
      continue;

    out << "cylinder {\n";
    out << std::format(
        "  <{}, {}, {}>, <{}, {}, {}>, CarbonBondCylinderRadius\n", atom1->x[0],
        atom1->x[1], atom1->x[2], atom2->x[0], atom2->x[1], atom2->x[2]);
    out << "  material { CarbonBondMaterial }\n";
    out << "}\n\n";
  }

  out << "// Carbon atoms\n";
  for (auto atom : get_atoms()) {
    if (atom->whether_located == false)
      continue;
    out << "sphere {\n";
    out << std::format("  <{}, {}, {}>, CarbonAtomSphereRadius\n", atom->x[0],
                       atom->x[1], atom->x[2]);
    out << "  material { CarbonAtomMaterial }\n";
    out << "}\n\n";
  }

  out << "// Carbon atoms\n";
  for (auto visit : get_vectors()) {
    auto a = visit->back_atom->x;
    auto b = visit->front_atom->x;
    auto capl = (b - a).normal() * 0.5;
    b -= capl;
    auto capend = b + capl * (1.0 - arrow_radius / 0.5);
    out << "merge {\n"
        << std::format(
               "        cylinder {{ <{}, {}, {}>, <{}, {}, {}>, "
               "CarbonArrowRadius}}\n",
               a[0], a[1], a[2], b[0], b[1], b[2])
        << std::format(
               "        cone {{<{}, {}, {}>, CarbonArrowRadius*2, <{}, {}, "
               "{}>, 0}}\n",
               b[0], b[1], b[2], capend[0], capend[1], capend[2])
        << std::format("        sphere {{<{}, {}, {}>, CarbonArrowRadius }}\n",
                       a[0], a[1], a[2]);
    if (visit->state == topo_state::i)
      out << "        material { ArrowStateI }\n";
    else if (visit->state == topo_state::j)
      out << "        material { ArrowStateJ }\n";
    else if (visit->state == topo_state::k)
      out << "        material { ArrowStateK }\n";
    else if (visit->state == topo_state::ni)
      out << "        material { ArrowStateNI }\n";
    else if (visit->state == topo_state::nj)
      out << "        material { ArrowStateNJ }\n";
    else if (visit->state == topo_state::nk)
      out << "        material { ArrowStateNK }\n";
    else
      out << "        material { CarbonAtomMaterial }\n";
    out << "}\n\n";
  }

  out.close();
}

template <typename T>
auto topology<T>::create_atom() -> atom_type* {
  auto atom = ::new (atom_pool.malloc()) atom_type();
  atoms.insert(atom);
  return atom;
}

template <typename T>
void topology<T>::delete_atom(atom_type* atom) {
  atoms.erase(atom);
  atom->~atom_type();
  atom_pool.free(atom);
}

template <typename T>
auto topology<T>::create_vector(topo_state state) -> topo_vector* {
  auto vector = ::new (vector_pool.malloc()) topo_vector(state);
  vectors.insert(vector);
  return vector;
}

template <typename T>
void topology<T>::delete_vector(topo_vector* vector) {
  vectors.erase(vector);
  vector->~topo_vector();
  vector_pool.free(vector);
}

template <typename T>
auto topology<T>::create_heptagon_path() -> topo_path_type {
  topo_path_type stack;
  stack.push_back(create_vector(topo_state::k));
  stack.push_back(create_vector(topo_state::ni));
  stack.push_back(create_vector(topo_state::nj));
  stack.push_back(create_vector(topo_state::nk));
  stack.push_back(create_vector(topo_state::i));
  stack.push_back(create_vector(topo_state::j));
  stack.push_back(create_vector(topo_state::k));
  for (auto i : number_range(7)) {
    stack[i]->link_front(stack[i + 1]);
  }
  for (auto i : number_range(7)) {
    stack[i]->set_front_atom(create_atom());
  }

  for (auto i : number_range(7)) {
    stack[i]->set_back_atom(stack[i - 1]->front_atom);
  }

  for (auto i : number_range(7)) {
    stack[i]->set_back_atom(stack[i - 1]->front_atom);
  }

  for (auto i : number_range(7)) {
    stack[i]->back_atom->connect_bond(stack[i]->front_atom);
  }

  DataTypeVector3D v(1.407 * std::sqrt(3.) / 2, 1.407 / 2, 0);
  for (auto i : number_range(7)) {
    stack[i]->back_atom->x = v;
    v.rotate({0, 0, 1}, M_PI / 3.5);
  }
  return stack;
}

template <typename T>
auto topology<T>::create_pentagon_path() -> topo_path_type {
  topo_path_type stack;
  stack.push_back(create_vector(topo_state::k));
  stack.push_back(create_vector(topo_state::ni));
  stack.push_back(create_vector(topo_state::nj));
  stack.push_back(create_vector(topo_state::nk));
  stack.push_back(create_vector(topo_state::i));
  for (auto i : number_range(0, 5)) {
    stack[i]->link_front(stack[i + 1]);
  }
  for (auto i : number_range(0, 5)) {
    stack[i]->set_front_atom(create_atom());
  }

  for (auto i : number_range(0, 5)) {
    stack[i]->set_back_atom(stack[i - 1]->front_atom);
  }

  for (auto i : number_range(0, 5)) {
    stack[i]->set_back_atom(stack[i - 1]->front_atom);
  }

  for (auto i : number_range(0, 5)) {
    stack[i]->back_atom->connect_bond(stack[i]->front_atom);
  }

  vector3<DataTypeFloat> v(1.307 * std::sqrt(3.) / 2, 1.307 / 2, 0);
  for (auto i : number_range(0, 5)) {
    stack[i]->back_atom->x = v;
    v.rotate({0, 0, 1}, M_PI / 2.5);
  }
  return stack;
}

template <typename T>
auto topology<T>::create_forward_path_to(const int n, const int m)
    -> topo_path_type {
  topo_path_type stack;

  if (n >= 0 and m >= 0) {
    std::stack<int> m_stack, n_stack;
    for (auto i : number_range(abs(m)))
      m_stack.push(1);
    for (auto i : number_range(abs(n)))
      n_stack.push(1);
    for (;;) {
      if (m_stack.size() == 0 and n_stack.size() == 0)
        break;
      if (n_stack.size()) {
        n_stack.pop();
        auto a = create_vector(topo_state::i);
        auto b = create_vector(topo_state::nk);
        a->link_front(b);
        b->link_back(a);
        if (stack.size()) {
          stack.back()->link_front(a);
        }
        stack.push_back(a);
        stack.push_back(b);
      }
      if (m_stack.size()) {
        m_stack.pop();
        auto a = create_vector(topo_state::i);
        auto b = create_vector(topo_state::j);
        a->link_front(b);
        b->link_back(a);
        if (stack.size()) {
          stack.back()->link_front(a);
        }
        stack.push_back(a);
        stack.push_back(b);
      }
    }
  } else if (n >= 1 and n + m >= 0) {
    std::stack<int> m_stack, n_stack;
    for (auto i : number_range(abs(m)))
      m_stack.push(1);
    for (auto i : number_range(abs(n) - abs(m)))
      n_stack.push(1);
    for (;;) {
      if (m_stack.size() == 0 and n_stack.size() == 0)
        break;
      if (n_stack.size()) {
        n_stack.pop();
        auto a = create_vector(topo_state::i);
        auto b = create_vector(topo_state::nk);
        a->link_front(b);
        b->link_back(a);
        if (stack.size()) {
          stack.back()->link_front(a);
        }
        stack.push_back(a);
        stack.push_back(b);
      }
      if (m_stack.size()) {
        int reverse = m_stack.top();
        m_stack.pop();
        auto a = create_vector(topo_state::nj);
        auto b = create_vector(topo_state::nk);
        a->link_front(b);
        b->link_back(a);
        if (stack.size()) {
          stack.back()->link_front(a);
        }
        stack.push_back(a);
        stack.push_back(b);
      }
    }
  } else
    error(np_debug_info, "Require n_i >= 0 and n_i+m_i >= 0!");

  DataTypeVector3D vi(1.407 * std::sqrt(3) / 2, 1.407 / 2, 0);
  DataTypeVector3D vj(0, 1.407, 0);
  DataTypeVector3D vnk(1.407 * std::sqrt(3) / 2, -1.407 / 2, 0);

  if (stack.size()) {
    stack[0]->set_back_atom(create_atom());
    stack[0]->back_atom->x = vi;
  }

  for (auto i : number_range(stack.size())) {
    stack[i]->set_front_atom(create_atom());
    if (stack[i]->front != nullptr) {
      stack[i]->front->set_back_atom(stack[i]->front_atom);
    }
    if (stack[i]->state == topo_state::i)
      stack[i]->front_atom->x = stack[i]->back_atom->x + vi;
    else if (stack[i]->state == topo_state::j)
      stack[i]->front_atom->x = stack[i]->back_atom->x + vj;
    else if (stack[i]->state == topo_state::nk)
      stack[i]->front_atom->x = stack[i]->back_atom->x + vnk;
    else if (stack[i]->state == topo_state::ni)
      stack[i]->front_atom->x = stack[i]->back_atom->x - vi;
    else if (stack[i]->state == topo_state::nj)
      stack[i]->front_atom->x = stack[i]->back_atom->x - vj;
    else if (stack[i]->state == topo_state::k)
      stack[i]->front_atom->x = stack[i]->back_atom->x - vnk;
    else
      error(np_debug_info, "Error happens when create a chain topo_vector");
  }
  for (auto i : number_range(int(stack.size()))) {
    stack[i]->back_atom->connect_bond(stack[i]->front_atom);
  }
  return stack;
}

template <typename T>
auto topology<T>::create_backward_path_of(topo_path_type& other)
    -> topo_path_type {
  topo_path_type stack;
  for (auto it = other.buff.rbegin(); it != other.buff.rend(); ++it) {
    auto tmp = create_vector((*it)->state);
    tmp->turn_back();
    tmp->set_front_atom((*it)->back_atom);
    tmp->set_back_atom((*it)->front_atom);
    stack.push_back(tmp);
  }
  const int size = (int)stack.size() - 1;
  for (int i = 0; i < size; i++) {
    stack[i]->link_front(stack[i + 1]);
  }
  return stack;
}

template <typename T>
auto topology<T>::remove_successive_reverse_vector(topo_vector*& start)
    -> bool {
  auto ret = false;
  for (auto current = start;;) {
    if (current->is_back(current->front)) {
      if (current->is_back_active() == false and
          current->front->is_front_active() == false) {
        std::cerr << "Overlap between rings!\n";
        std::cerr
            << "!  Pentagonal carbon rings are not allowed to overlap with "
               "each other.\n"
               "!  This cased by placing two pentagon rings at same position!\n"
               "!  e.g., [..(0,1),(1,0)..] or [(0,0),(1,0)...(0,1)]!\n";
        error(np_debug_info, "Runtime error!");
      }
      auto pax = (current->front->front_atom->x + current->back_atom->x) / 2;
      start = current->back;
      start->link_front(current->front->front);
      auto atom = start->front_atom;
      start->front_atom->self_replaced_by(start->front->back_atom);
      delete_atom(atom);
      start->front_atom->x = pax;
      delete_vector(current->front);
      delete_vector(current);
      current = start;
      ret = true;
      continue;
    }
    current = current->front;
    if (current == start)
      break;
  }
  return ret;
}

template <typename T>
auto topology<T>::get_number_of_vectors_between(topo_vector* a,
                                                topo_vector* b) const
    -> size_t {
  size_t count = 0;
  for (auto current = a;;) {
    if (current == b)
      break;
    current = current->front;
    if (current == nullptr)
      error(np_debug_info,
            "Unexpected break in topological path while calculating vector "
            "distance.\n"
            "Graphene-derived nanostructures present significant computational "
            "challenges.");
    count++;
  }
  return count;
}

template <typename T>
auto topology<T>::find_most_active_site(topo_vector* start) const
    -> topo_path_type {
  topo_path_type ret;
  size_t ret_num = 0;
  start = start->find_next_front_activity();

  for (auto current = start;;) {
    circular_topo_path<topo_vector*> result;
    size_t current_number = 1;
    result.push_back(current);
    auto next = current->front;
    for (;;) {
      result.push_back(next);
      if (next->is_front_active())
        break;
      next = next->front;
      current_number++;
    }
    if (current_number > ret_num) {
      ret = result;
      ret_num = current_number;
    } else if (current_number == ret_num) {
      if (result[0]->front_atom->x.length() < ret[0]->front_atom->x.length())
        ret = result;
    }
    current = next;
    if (current == start)
      break;
  }
  return ret;
}

template <typename T>
auto topology<T>::find_active_site_contains(topo_vector* st) const
    -> topo_path_type {
  topo_path_type result;
  auto tmpa = st->find_previous_front_activity();
  auto tmpb = st->find_next_front_activity();
  result.push_back(tmpa);
  result.push_back(tmpb);
  return result;
}

template <typename T>
auto topology<T>::get_number_on_the_edgepath(topo_vector* start) const
    -> size_t {
  size_t count = 0;
  auto current = start;
  for (;;) {
    current = current->front;
    count++;
    if (current == start)
      break;
  }
  return count;
}

template <typename T>
auto topology<T>::add_hexagon_onto_active_site(topo_path_type& path)
    -> topo_vector* {
  const auto path_len = path.size();
  switch (path_len) {
    case 3: {
      auto pax = path[0]->front_atom->third_nei_coordinate();
      auto pcx = path[-1]->front_atom->third_nei_coordinate();
      auto tempa = path[1];
      tempa->set_state(path[0]->state);
      tempa->turn_right();
      auto tempb = create_vector(tempa->state);
      tempb->turn_left();
      tempa->link_front(tempb);
      auto tempc = create_vector(tempb->state);
      tempc->turn_left();
      tempb->link_front(tempc);
      auto tempd = path[-1];
      tempd->set_state(tempc->state);
      tempd->turn_left();
      tempc->link_front(tempd);

      tempa->set_front_atom(create_atom());
      tempb->set_back_atom(tempa->front_atom);

      tempb->set_front_atom(create_atom());
      tempc->set_back_atom(tempb->front_atom);

      tempc->set_front_atom(create_atom());
      tempd->set_back_atom(tempc->front_atom);

      tempa->front_atom->connect_bond(tempa->back_atom);
      tempb->front_atom->connect_bond(tempb->back_atom);
      tempc->front_atom->connect_bond(tempc->back_atom);
      tempd->front_atom->connect_bond(tempd->back_atom);

      tempa->front_atom->x = pax;
      tempc->front_atom->x = pcx;
      tempb->front_atom->x = (pax + pcx) / 2 +
                             ((tempa->front_atom->x - tempa->back_atom->x) -
                              (tempd->front_atom->x - tempd->back_atom->x)).normal()*0.7;

      return tempa;
    }
    case 4: {
      auto pax = path[0]->front_atom->third_nei_coordinate();
      auto pcx = path[-1]->front_atom->third_nei_coordinate();
      auto tempa = path[1];
      tempa->set_state(path[0]->state);
      tempa->turn_right();
      auto tempb = path[2];
      tempb->set_state(tempa->state);
      tempb->turn_left();
      auto tempc = path[3];
      tempc->set_state(tempb->state);
      tempc->turn_left();

      tempa->set_front_atom(create_atom());
      tempb->set_back_atom(tempa->front_atom);

      tempb->set_front_atom(create_atom());
      tempc->set_back_atom(tempb->front_atom);

      tempa->front_atom->connect_bond(tempa->back_atom);
      tempb->front_atom->connect_bond(tempb->back_atom);
      tempc->front_atom->connect_bond(tempc->back_atom);

      tempa->front_atom->x = pax;
      tempb->front_atom->x = pcx;

      return tempa;
    }
    case 5: {
      auto pax = path[0]->front_atom->third_nei_coordinate();
      auto pcx = path[-1]->front_atom->third_nei_coordinate();
      auto tempa = path[1];
      tempa->set_state(path[0]->state);
      tempa->turn_right();
      auto tempb = path[-1];
      tempb->set_state(tempa->state);
      tempb->turn_left();
      tempa->link_front(tempb);
      tempa->set_front_atom(create_atom());
      tempb->set_back_atom(tempa->front_atom);
      tempa->front_atom->connect_bond(tempa->back_atom);
      tempb->front_atom->connect_bond(tempb->back_atom);
      delete_vector(path[2]);
      delete_vector(path[3]);

      tempa->front_atom->x = (pax + pcx) / 2;
      return tempa;
    }
    case 6: {
      auto pax = path[0]->front_atom->third_nei_coordinate();
      auto pcx = path[-1]->front_atom->third_nei_coordinate();
      auto tempa = path[1];
      tempa->set_state(path[0]->state);
      tempa->turn_right();
      tempa->link_front(path[-1]->front);
      tempa->set_front_atom(path[-1]->front_atom);
      tempa->front_atom->connect_bond(tempa->back_atom);
      delete_vector(path[2]);
      delete_vector(path[3]);
      delete_vector(path[4]);
      delete_vector(path[5]);
      return tempa;
    }
    case 7: {
      auto atom = path[0]->front_atom;
      atom->self_replaced_by(path[-1]->front_atom);
      path[0]->link_front(path[-1]->front);
      delete_atom(atom);
      delete_vector(path[1]);
      delete_vector(path[2]);
      delete_vector(path[3]);
      delete_vector(path[4]);
      delete_vector(path[5]);
      delete_vector(path[6]);
      return path[0]->front;
    }
    default: {
      std::cerr << "Overlap between rings!\n";
      std::cerr
          << "!  Pentagonal carbon rings are not allowed to overlap with each "
             "other.\n"
             "!  This cased by placing two pentagon rings at same position!\n"
             "!  e.g., [..(0,1),(1,0)..] or [(0,0),(1,0)...(0,1)]!\n";
      error(np_debug_info,
            std::format(
                "Overlap error: detected {} edges between two neighboring "
                "active atoms,\nwhich may be caused by overlapping between "
                "pentagonal carbon rings during the earlier process\n",
                path_len - 1));
      return nullptr;  // less warnings
    }
  }
}

template <typename T>
auto topology<T>::dismantle_hexagons_from_edge_path(topo_vector* _start)
    -> topo_vector* {
  auto start = _start->find_next_front_activity();

  for (auto current = start;;) {
    current = current->find_next_front_activity();
    auto next = current->front->find_next_front_activity();

    if (next == current->front) {
      auto v1 = current->back;
      auto v2 = v1->front;
      auto v3 = v2->front;
      auto v4 = v3->front;
      auto v5 = v4->front;

      auto atom_a = v1->front_atom;
      auto atom_b = v4->front_atom;

      atom_type* atom_c = nullptr;
      atom_type* atom_d = nullptr;
      for (auto t : atom_a->conn_vertex) {
        if (t != v1->back_atom and t != v2->front_atom) {
          atom_c = t;
          break;
        }
      }
      for (auto t : atom_b->conn_vertex) {
        if (t != v4->back_atom and t != v5->front_atom) {
          atom_d = t;
          break;
        }
      }

      if (atom_c != nullptr and atom_c->is_connected(atom_d)) {
        delete_atom(v2->front_atom);
        delete_atom(v3->front_atom);
        v2->set_front_atom(atom_c);
        v3->set_back_atom(atom_c);
        v3->set_front_atom(atom_d);
        v4->set_back_atom(atom_d);
        v2->state = v1->state;
        v2->turn_left();
        v3->state = v2->state;
        v3->turn_right();
        v4->state = v3->state;
        v4->turn_right();
        start = current->find_next_front_activity();
        current = start;
        continue;
      }
    }
    current = next;
    if (current == start)
      break;
  }
  return start;
}

template <typename T>
auto topology<T>::get_sum_from_path(topo_vector* start) const
    -> std::array<int, 3> {
  auto current = start;
  int i{0}, j{0}, k{0};
  for (;;) {
    switch (current->state) {
      case topo_state::i:
        i++;
        break;
      case topo_state::j:
        j++;
        break;
      case topo_state::k:
        k++;
        break;
      case topo_state::ni:
        i--;
        break;
      case topo_state::nj:
        j--;
        break;
      case topo_state::nk:
        k--;
        break;
    }
    current = current->front;
    if (current == start or current == nullptr)
      break;
  }

  return {i, j, k};
}

}  // namespace grapology

#endif  // GRAPOLOGY_TOPOLOGY_HPP