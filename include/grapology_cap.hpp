// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGYCARBONPOLOGY
#define GRAPOLOGYCARBONPOLOGY

namespace grapology {
/**
 * @class cap
 * @brief Class for creating and manipulating carbon nanotubes caps
 *
 * This class provides the fundamental functionality for generating carbon
 * nanotube caps.
 */

class carbon_cap : public topology<carbon_cap> {
 public:
  bool whether_optimizing;
  explicit carbon_cap(bool whether_optimizing = true)
      : topology<carbon_cap>(), whether_optimizing(whether_optimizing) {}
  ~carbon_cap() {}
  /**
   * @brief Constructs a derived structure based on the given parameters
   *
   * This function serves as the implementation of the construct() function in
   * the basis topology class. It sequentially performs two key steps:
   * 1. Creates the carbon nanotube cap structure
   * 2. Dynamically grows from the cap edge to form the carbon nanotube body
   *
   * The function receives all necessary parameters through an integer vector,
   * which define the chirality and growth characteristics of the carbon
   * nanotube.
   *
   * @param mn Integer vector containing construction parameters, must include
   * at least 13 elements:
   *           - mn[0]-mn[11]: First to 12th chiral indices n1,m1-n6,m6
   *           - mn[12]: Number of hexagons to add, controlling the length of
   * the carbon nanotube
   */
  void construct_derived_structure(const std::vector<int>& mn) {
    if (mn.size() < 13) {
      error(np_debug_info, "Input vector must contain at least 13 elements");
    }
    auto visit = create_cap(mn[0], mn[2], mn[4], mn[6], mn[8], mn[10], 
                            mn[1], mn[3], mn[5], mn[7], mn[9], mn[11]);
    visit = dynamic_growth(visit, mn[12]);
  }

 private:
  /**
   * @brief Creates a carbon nanotube cap based on specified chiral indices of
   * the cap
   *
   * This method generates a carbon nanotube cap structure using the provided
   * chiral indices which define the geometry and chirality of the cap. It
   * constructs the cap by creating and connecting topological vectors that form
   * the hexagonal and pentagonal rings characteristic of carbon nanotube caps.
   * The method ensures that the generated cap is structurally stable and
   * adheres to the constraints of the carbon nanotube cap geometry.
   *
   * @param n1-n6, m1-m6 Chiral indices defining the cap structure. These
   * parameters control the geometry, size, and chirality of the resulting
   * carbon nanotube cap. The sum of n1-n6 and m1-m6 determine the chirality of
   * the later grown tube.
   *
   * @return Pointer to a topological vector on the edge of the cap structure,
   * which can be used as a starting point for further growth of the carbon
   * nanotube.
   *
   * @note The function performs structural optimization if the
   * whether_optimizing flag is set to true (default). This involves dynamic
   * relaxation steps to ensure the resulting structure is energetically stable.
   *
   * @warning If two pentagons are placed at same position, the function will
   * throw an error. Similarly, if an error occurs during the folding process
   * due to ring overlapping, an error will be raised.
   *
   * The creation process involves several key steps:
   * 1. Creating vector chains (path) based on the provided chiral indices
   * 2. Connecting these chains (path) to form a closed path
   * 3. Removing overlaps and optimizing the structure
   * 4. Adding hexagons to complete the cap geometry
   * 5. Final optimization and spherical position folding
   */
  auto create_cap(int n1, int n2, int n3, int n4, int n5, int n6,
                  int m1, int m2, int m3, int m4, int m5, int m6) -> topo_vector* {
    const int cn = n1 + n2 + n3 + n4 + n5 + n6;
    const int cm = m1 + m2 + m3 + m4 + m5 + m6;

    auto create_vector_based_chain =
        [&](int n, int m,
            DataTypeFloat angle) -> circular_topo_path<topo_vector*> {
      auto pentagon = create_pentagon_path();
      pentagon.move_atoms(n, m);
      pentagon.rotate_atoms(angle);

      auto forward_v = create_forward_path_to(n, m);
      forward_v.rotate_atoms(angle);
      auto backward_v = create_backward_path_of(forward_v);
      backward_v.turn_vector_right(1);
      auto tmp = create_vector(topo_state::j);
      if (backward_v.size())
        backward_v[-1]->link_front(tmp);
      if (backward_v.size())
        tmp->set_back_atom(backward_v[-1]->front_atom);
      else
        tmp->set_back_atom(pentagon[-1]->front_atom);
      backward_v.push_back(tmp);
      if (forward_v.size()) {
        auto atom = forward_v[-1]->front_atom;
        atom->self_replaced_by(pentagon[0]->back_atom);
        delete_atom(atom);
      }
      return forward_v + pentagon + backward_v;
    };

    std::vector<circular_topo_path<topo_vector*>> chains_v;
    chains_v.reserve(6);

    DataTypeFloat angle =
        (n1 | m1 and n2 | m2 and n3 | m3 and n4 | m4 and n5 | m5 and n6 | m6)
            ? static_cast<DataTypeFloat>(M_PI / 3.0)
            : static_cast<DataTypeFloat>(M_PI / 2.5);

    int intial_rotate_angle_index = 0;

    if (n1 != 0 or m1 != 0) {
      chains_v.push_back(
          create_vector_based_chain(n1, m1, intial_rotate_angle_index * angle));
      intial_rotate_angle_index++;
    }
    if (n2 != 0 or m2 != 0) {
      chains_v.push_back(
          create_vector_based_chain(n2, m2, intial_rotate_angle_index * angle));
      intial_rotate_angle_index++;
    }
    if (n3 != 0 or m3 != 0) {
      chains_v.push_back(
          create_vector_based_chain(n3, m3, intial_rotate_angle_index * angle));
      intial_rotate_angle_index++;
    }
    if (n4 != 0 or m4 != 0) {
      chains_v.push_back(
          create_vector_based_chain(n4, m4, intial_rotate_angle_index * angle));
      intial_rotate_angle_index++;
    }
    if (n5 != 0 or m5 != 0) {
      chains_v.push_back(
          create_vector_based_chain(n5, m5, intial_rotate_angle_index * angle));
      intial_rotate_angle_index++;
    }
    if (n6 != 0 or m6 != 0) {
      chains_v.push_back(
          create_vector_based_chain(n6, m6, intial_rotate_angle_index * angle));
      intial_rotate_angle_index++;
    }

    for (size_t i = 0; i < chains_v.size() - 1; i++) {
      chains_v[i][-1]->set_front_atom(chains_v[i + 1][0]->back_atom);
      chains_v[i][-1]->front_atom->connect_bond(chains_v[i][-1]->back_atom);
    }
    chains_v.back()[-1]->set_front_atom(chains_v.front()[0]->back_atom);
    chains_v.back()[-1]->front_atom->connect_bond(
        chains_v.back()[-1]->back_atom);

    if (chains_v.size() == 6) {
      auto tmp = chains_v[0] + chains_v[1] + chains_v[2] + chains_v[3] +
                 chains_v[4] + chains_v[5] + chains_v[0];
    } else if (chains_v.size() == 5) {
      auto tmp = chains_v[0] + chains_v[1] + chains_v[2] + chains_v[3] +
                 chains_v[4] + chains_v[0];
    } else
      error(np_debug_info, "Two pentagons cannot be at (0,0)!");

    auto visit = chains_v[0][0];

    locate loc(get_atoms());
    if (remove_successive_reverse_vector(visit) and whether_optimizing)
      loc.process_relaxing();

    for (;;) {
      auto&& actives = find_most_active_site(visit);
      const auto num_vectors = actives.size();
      if (num_vectors <= 4)
        break;
      visit = add_hexagon_onto_active_site(actives);
      if (whether_optimizing)
        loc.process_relaxing();
      if (remove_successive_reverse_vector(visit) and whether_optimizing)
        loc.process_relaxing();
    }

    //if (remove_successive_reverse_vector(visit) and whether_optimizing)
    //  loc.process_relaxing();

    loc.spherical_position_folding();

    for (const size_t bonds_on_chain_st = get_number_on_the_edgepath(visit);;) {
      const size_t bonds_on_chain = get_number_on_the_edgepath(visit);
      if (bonds_on_chain <= 2 * (cn + cm))
        break;
      else if (bonds_on_chain > bonds_on_chain_st) {
        std::cerr << "Overlap between rings!\n";
        std::cerr
            << "!  Pentagonal carbon ring is not allowed to overlap with each "
               "other.\n"
               "!  This cased by placing two pentagon rings at same position!\n"
               "!  e.g., [..(0,1),(1,0)..] or [(0,0),(1,0)...(0,1)]!\n";
        error(np_debug_info, "Runtime error!");
      }
      auto&& actives = find_most_active_site(visit);
      visit = add_hexagon_onto_active_site(actives);
      if (whether_optimizing)
        loc.process_relaxing();
      if (remove_successive_reverse_vector(visit) and whether_optimizing)
        loc.process_relaxing();
    }

    return dismantle_hexagons_from_edge_path(visit);
  }

  /**
   * @brief Grows dynamically from a carbon cap edge to a carbon nanotube
   * @param start A starting topological vector (on the edge of the cap) for
   * growth
   * @param hexagon_num Target size, amount of hexagons to add to the grown
   * structure, which controls the length of the grown structure.
   * @return Pointer to the final topological vector path of the grown
   * structure, which can be used for analysis or further growth operations.
   */
  auto dynamic_growth(topo_vector* start, const size_t hexagon_num)
      -> topo_vector* {
    size_t count_hexagon{0};
    locate loc(get_atoms());
    locate<false> loc2(get_atoms());

    auto current = start;
    for (;;) {
      if (++count_hexagon > hexagon_num)
        break;
      auto&& actives = find_most_active_site(current);
      current = add_hexagon_onto_active_site(actives);
      if (whether_optimizing)
        loc.process_relaxing();
      if (remove_successive_reverse_vector(current) and whether_optimizing)
        loc.process_relaxing();
    }

    if (whether_optimizing)
      loc2.dynamic_relaxing(100);
    return current;
  }

  /**
   * @brief Calculates the chiral indices of a carbon nanotube from an edge path
   *
   * This function computes the chiral indices of a carbon nanotube, which are
   * important parameters that characterize the structure of carbon nanotubes.
   * Chiral indices are typically represented as (n,m) and determine the
   * electrical and mechanical properties of the nanotube. The function extracts
   * and calculates these indices by analyzing the topological vectors of a
   * given edge path.
   *
   * @param start Pointer to the starting topological vector of the edge path
   * @return An array containing two integers representing the computed chiral
   * indices (n,m) The returned chiral indices ensure n ≥ m, which is the
   * standard representation for carbon nanotube chiral indices
   *
   * @note Chiral indices determine the type of carbon nanotube:
   *       - When m = 0, it's a zigzag carbon nanotube
   *       - When n = m, it's an armchair carbon nanotube
   *       - In other cases, it's a chiral carbon nanotube
   */
  auto cal_chiral_indices_from_edgepath(topo_vector* start)
      const -> std::array<int, 2>  // return the chiral index for checking the path
  {
    auto sum = get_sum_from_path(start);
    const int i = sum[0];
    const int j = sum[1];
    const int k = sum[2];
    if (j != i + k)
      error(np_debug_info, "Chiral indices are not valid!");
    return {i, k};
  }

}; // class carbon_cap

}  // namespace grapology

#endif
