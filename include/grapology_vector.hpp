// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGYVECTOR
#define GRAPOLOGYVECTOR

namespace grapology {

/**
 * @enum topo_state
 * @brief Enumerates the possible states of a topological vector
 *
 * These states represent different positions in the hexagonal lattice structure
 * of graphene. Each state corresponds to a relative orientation or direction
 * within the lattice.

             k
          <_____
       ni/      \j
        /        \
      nj\        /i
         \_____>/
            nk
*/
enum class topo_state { i, j, k, ni, nj, nk };  // enum class

/**
 * @class topo_vector
 * @brief Represents a topological vector in graphene structures
 *
 * This class encapsulates the topological information of vectors connecting
 * atoms in graphene-based structures. It manages vector orientations,
 * connections to adjacent vectors, and associations with atoms.
 * The class provides methods for manipulating vector directions and
 * establishing connections in the lattice structure.
 */

struct alignas(8) topo_vector final {
  using atom_type = vertex<topo_vector>;
  topo_vector *front, *back;
  atom_type *front_atom, *back_atom;
  topo_state state;

  inline void link_front(topo_vector* f) {
    if (front != f) {
      break_front();
      front = f;
      f->break_back();
      f->back = this;
    }
  }
  inline void link_back(topo_vector* b) {
    if (back != b) {
      break_back();
      back = b;
      b->break_front();
      b->front = this;
    }
  }
  inline void break_front() {
    if (front) {
      if (front->back == this)
        front->back = nullptr;
      front = nullptr;
    }
  }
  inline void break_back() {
    if (back) {
      if (back->front == this)
        back->front = nullptr;
      back = nullptr;
    }
  }

  inline void turn_left() {
    switch (state) {
      case topo_state::i:
        state = topo_state::j;
        break;
      case topo_state::j:
        state = topo_state::k;
        break;
      case topo_state::k:
        state = topo_state::ni;
        break;
      case topo_state::ni:
        state = topo_state::nj;
        break;
      case topo_state::nj:
        state = topo_state::nk;
        break;
      default:
        state = topo_state::i;
        break;
    }
  }

  inline void turn_left(size_t n) {
    for (auto i = 0; i < n; ++i)
      turn_left();
  }

  inline void turn_right() {
    switch (state) {
      case topo_state::i:
        state = topo_state::nk;
        break;
      case topo_state::j:
        state = topo_state::i;
        break;
      case topo_state::k:
        state = topo_state::j;
        break;
      case topo_state::ni:
        state = topo_state::k;
        break;
      case topo_state::nj:
        state = topo_state::ni;
        break;
      default:
        state = topo_state::nj;
        break;
    }
  }

  inline void turn_right(size_t n) {
    for (auto i = 0; i < n; ++i)
      turn_right();
  }

  inline void turn_back() {
    switch (state) {
      case topo_state::i:
        state = topo_state::ni;
        break;
      case topo_state::j:
        state = topo_state::nj;
        break;
      case topo_state::k:
        state = topo_state::nk;
        break;
      case topo_state::ni:
        state = topo_state::i;
        break;
      case topo_state::nj:
        state = topo_state::j;
        break;
      default:
        state = topo_state::k;
        break;
    }
  }

  inline void turn_back(size_t n) {
    for (auto i = 0; i < n; ++i)
      turn_back();
  }

  inline auto is_right(const topo_vector* other) const -> bool {
    switch (state) {
      case topo_state::i:
        if (other->state == topo_state::nk)
          return true;
        break;
      case topo_state::j:
        if (other->state == topo_state::i)
          return true;
        break;
      case topo_state::k:
        if (other->state == topo_state::j)
          return true;
        break;
      case topo_state::ni:
        if (other->state == topo_state::k)
          return true;
        break;
      case topo_state::nj:
        if (other->state == topo_state::ni)
          return true;
        break;
      default:
        if (other->state == topo_state::nj)
          return true;
        break;
    }
    return false;
  }

  inline auto is_left(const topo_vector* other) const -> bool {
    switch (state) {
      case topo_state::i:
        if (other->state == topo_state::j)
          return true;
        break;
      case topo_state::j:
        if (other->state == topo_state::k)
          return true;
        break;
      case topo_state::k:
        if (other->state == topo_state::ni)
          return true;
        break;
      case topo_state::ni:
        if (other->state == topo_state::nj)
          return true;
        break;
      case topo_state::nj:
        if (other->state == topo_state::nk)
          return true;
        break;
      default:
        if (other->state == topo_state::i)
          return true;
        break;
    }
    return false;
  }

  inline auto is_back(const topo_vector* other) const -> bool {
    switch (state) {
      case topo_state::i:
        if (other->state == topo_state::ni)
          return true;
        break;
      case topo_state::j:
        if (other->state == topo_state::nj)
          return true;
        break;
      case topo_state::k:
        if (other->state == topo_state::nk)
          return true;
        break;
      case topo_state::ni:
        if (other->state == topo_state::i)
          return true;
        break;
      case topo_state::nj:
        if (other->state == topo_state::j)
          return true;
        break;
      default:
        if (other->state == topo_state::k)
          return true;
        break;
    }
    return false;
  }

  inline auto set_front_atom(atom_type* atom) -> bool {
    if (front_atom == atom)
      return false;
    if (front_atom) {
      front_atom->conn_vector.erase(this);
    }

    front_atom = atom;
    front_atom->conn_vector.insert(this);

    return true;
  }

  inline auto set_back_atom(atom_type* atom) -> bool {
    if (back_atom == atom)
      return false;
    if (back_atom) {
      back_atom->conn_vector.erase(this);
    }

    back_atom = atom;
    back_atom->conn_vector.insert(this);

    return true;
  }

  inline auto is_front_active() const -> bool { return is_left(front); }

  inline auto is_back_active() const -> bool { return back->is_left(this); }

  inline auto find_next_front_activity() -> topo_vector* {
    auto current = this;
    while (true) {
      if (current->is_front_active())
        return current;
      current = current->front;
      if (current == this)
        error(np_debug_info, "No active position found");
    }
    error(np_debug_info, "No active position found");
    return current;  // less warnings
  }

  inline auto find_previous_front_activity() -> topo_vector* {
    auto current = this;
    while (true) {
      if (current->is_front_active())
        return current;
      current = current->back;
      if (current == this)
        error(np_debug_info, "No active position found");
    }
    error(np_debug_info, "No active position found");
    return current;  // less warnings
  }

  void set_state(topo_state b) { state = b; }

  explicit topo_vector(topo_state b)
      : front(nullptr),
        back(nullptr),
        front_atom(nullptr),
        back_atom(nullptr),
        state(b) {}
  ~topo_vector() {
    if (front_atom)
      front_atom->conn_vector.erase(this);
    if (back_atom)
      back_atom->conn_vector.erase(this);
    if (front and front->back == this)
      front->back = nullptr;
    if (back and back->front == this)
      back->front = nullptr;
  }
  topo_vector(topo_vector&) = delete;
  topo_vector(topo_vector&&) = delete;
  static void * operator new(size_t) = delete;
  static void operator delete(void*) = delete;
};  // class topologic vector base

/**
 * @class circular_topo_path
 * @brief Template class for managing circular paths of topological vectors
 *
 * This class provides a container for topological vectors that form a circular
 * (closed-loop) path. It offers various operations for manipulating the path,
 * including vector rotation, atom movement, and path concatenation. The class
 * implements circular indexing semantics, allowing wrap-around access to
 * elements.
 *
 * @tparam T Type of elements stored in the path (defaults to topo_vector*)
 */

template <typename T = topo_vector*>
struct circular_topo_path {
  std::vector<T> buff;

  circular_topo_path() {}
  circular_topo_path(const circular_topo_path& other) : buff(other.buff) {}

  template <typename INT>
  auto operator[](const INT ind) -> T& {
    if constexpr (std::is_same_v<INT, int>) {
      if (ind < 0) {
        int size = (int)buff.size();
        return buff[ind % size == 0 ? 0 : ind % size + size];
      } else if (ind >= buff.size())
        return buff[ind % buff.size()];
      else
        return buff[ind];
    } else if constexpr (std::is_same_v<INT, size_t>) {
      if (ind >= buff.size())
        return buff[ind % buff.size()];
      else
        return buff[ind];
    } else {
      static_assert(false, "invalid index type : index must be int or size_t");
    }
  }

  void push_back(T v) { buff.push_back(v); }

  void pop_back() { buff.pop_back(); }

  auto front() -> T { return buff.front(); }

  auto back() -> T { return buff.back(); }

  auto size() const -> size_t { return buff.size(); }

  void turn_vector_left(size_t n = 1) {
    for (auto i : buff) {
      i->turn_left(n);
    }
  }

  void turn_vector_right(size_t n = 1) {
    for (auto i : buff) {
      i->turn_right(n);
    }
  }

  void rotate_atoms(DataTypeFloat angle) {
    std::set<vertex<topo_vector>*> con;
    for (auto i : buff) {
      if (i->front_atom)
        con.insert(i->front_atom);
      if (i->back_atom)
        con.insert(i->back_atom);
    }
    for (auto i : con)
      i->x.rotate({0, 0, 1}, angle);
  }

  void move_atoms(int n, int m) {
    std::set<vertex<topo_vector>*> con;
    for (auto i : buff) {
      if (i->front_atom)
        con.insert(i->front_atom);
      if (i->back_atom)
        con.insert(i->back_atom);
    }
    DataTypeVector3D tmp =
        DataTypeVector3D{static_cast<DataTypeFloat>(1.407 * std::sqrt(3.)), 0, 0} * n +
        DataTypeVector3D{static_cast<DataTypeFloat>(1.407 * std::sqrt(3.) / 2), static_cast<DataTypeFloat>(1.407 * 1.5), 0} * m;
    for (auto i : con)
      i->x += tmp;
  }

  auto operator+(circular_topo_path& other) -> circular_topo_path {
    circular_topo_path ret;
    for (auto i : this->buff)
      ret.buff.push_back(i);
    for (auto i : other.buff)
      ret.buff.push_back(i);
    if (buff.size() and other.buff.size())
      this->buff.back()->link_front(other.buff.front());
    return ret;
  }
};  // end of class topologic vertor container

template <typename INT>
struct number_range  // replace of std::range::iota from c++23
{
  using value_type = INT;
  struct GForwardIteratorEnd {};

  struct GForwardIterator {
    using iterator_category = std::forward_iterator_tag;
    value_type start, ends, increment, current;
    explicit constexpr GForwardIterator(value_type st,
                                        value_type en,
                                        value_type inc)
        : start(st), ends(en), increment(inc), current(st) {
      if (inc == 1 and en < st)
        increment = -1;
    }

    constexpr inline auto operator*() const -> value_type const  { return current; }
    constexpr inline auto operator++() -> GForwardIterator& {
      ++current;
      return *this;
    }
    constexpr inline auto operator++(int) -> GForwardIterator& {
      ++current;
      return *this;
    }
    friend constexpr inline bool operator==(const GForwardIterator& my,
                                            const GForwardIteratorEnd& other) {
      if (my.current == my.ends)
        return true;
      else
        return false;
    }
    friend constexpr inline bool operator!=(const GForwardIterator& my,
                                            const GForwardIteratorEnd& other) {
      return not(my == other);
    }
  };

  INT start, ends, increment;
  constexpr number_range(INT st, INT en, INT inc)
      : start(st), ends(en), increment(inc) {}
  constexpr number_range(INT st, INT en) : start(st), ends(en), increment(1) {}
  constexpr number_range(INT en) : start(0), ends(en), increment(1) {}

  constexpr inline auto begin() -> GForwardIterator { return GForwardIterator(start, ends, increment); }
  constexpr inline auto end() -> GForwardIteratorEnd { return GForwardIteratorEnd(); }
};  // end of class

}  // namespace grapology

#endif
