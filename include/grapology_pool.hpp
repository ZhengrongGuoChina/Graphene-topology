// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGY_POOL_HPP
#define GRAPOLOGY_POOL_HPP

namespace grapology {

/**
 * Fixed-size memory pool implementation for vertex and topo_vector
 * Allocates fixed-size memory blocks for cache efficiency
 */
struct grapology_pool {
 private:
  struct MemoryBlock {
    MemoryBlock* next;
  };

  const size_t block_size_;
  const size_t chunk_size_;
  MemoryBlock* free_list_;
  std::vector<void*> chunks_;

 public:
  explicit grapology_pool(size_t block_size, size_t chunk_size)
      : block_size_(std::max(sizeof(MemoryBlock), block_size)),
        chunk_size_(chunk_size),
        free_list_(nullptr) {}

  ~grapology_pool() {
    for (void* chunk : chunks_) {
      free(chunk);
    }
  }

  void* malloc() {
    if (free_list_ == nullptr) {
      void* chunk = new char[block_size_ * chunk_size_];
      if (chunk == nullptr) {
        error(np_debug_info,
              "grapology_pool::malloc: failed to allocate memory");
        return nullptr;
      }

      chunks_.push_back(chunk);

      char* block = static_cast<char*>(chunk);
      for (size_t i = 0; i < chunk_size_ - 1; ++i) {
        MemoryBlock* curr =
            reinterpret_cast<MemoryBlock*>(block + i * block_size_);
        curr->next =
            reinterpret_cast<MemoryBlock*>(block + (i + 1) * block_size_);
      }

      MemoryBlock* last = reinterpret_cast<MemoryBlock*>(
          block + (chunk_size_ - 1) * block_size_);
      last->next = nullptr;

      free_list_ = reinterpret_cast<MemoryBlock*>(block);
    }

    MemoryBlock* block = free_list_;
    free_list_ = block->next;
    return block;
  }

  void free(void* ptr) {
    if (ptr == nullptr) {
      return;
    }

    MemoryBlock* block = static_cast<MemoryBlock*>(ptr);
    block->next = free_list_;
    free_list_ = block;
  }

  grapology_pool(const grapology_pool&) = delete;
  grapology_pool& operator=(const grapology_pool&) = delete;
  grapology_pool(grapology_pool&&) = delete;
  grapology_pool& operator=(grapology_pool&&) = delete;
};

}  // namespace grapology

#endif  // GRAPOLOGY_POOL_HPP