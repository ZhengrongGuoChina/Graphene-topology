#ifndef GRAPOLOGY_SMALL_SET_HPP
#define GRAPOLOGY_SMALL_SET_HPP

#include <algorithm>
#include <cstddef>
#include <iterator>

namespace grapology {

/**
 * Unordered set that can store up to 6 elements, with extreme optimizations for small-scale data
 * Compatible with the basic interface of std::unordered_set, but avoids the overhead of hash tables
 */
template <typename T, size_t MAX_SIZE = 6>
struct small_unordered_set {
private:
    static_assert(std::is_copy_constructible_v<T>, "Type must be copy constructible");
    static_assert(std::is_copy_assignable_v<T>, "Type must be copy assignable");
    
    T data_[MAX_SIZE];
    
    size_t size_ = 0;
    
    constexpr T* find_element(const T& value) noexcept {
        for (size_t i = 0; i < size_; ++i) {
            if (data_[i] == value) {
                return &data_[i];
            }
        }
        return nullptr;
    }
    
    constexpr const T* find_element(const T& value) const noexcept {
        for (size_t i = 0; i < size_; ++i) {
            if (data_[i] == value) {
                return &data_[i];
            }
        }
        return nullptr;
    }

public:
    class iterator {
    private:
        T* ptr_;
        small_unordered_set& container_;
        size_t index_;
        
        friend class small_unordered_set;
        
        constexpr iterator(T* ptr, small_unordered_set& container, size_t index) noexcept
            : ptr_(ptr), container_(container), index_(index) {}

    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using pointer = T*;
        using reference = T&;
        
        constexpr iterator() = delete;
        
        constexpr reference operator*() const noexcept { return *ptr_; }
        constexpr pointer operator->() const noexcept { return ptr_; }
        
        constexpr iterator& operator++() noexcept {
            if (index_ + 1 < container_.size_) {
                ptr_ = &(container_.data_[++index_]);
            } else {
                ptr_ = nullptr;
                index_ = container_.size_;
            }
            return *this;
        }
        
        constexpr iterator& operator++(int) noexcept {
            if (index_ + 1 < container_.size_) {
                ptr_ = &(container_.data_[++index_]);
            } else {
                ptr_ = nullptr;
                index_ = container_.size_;
            }
            return *this;
        }
        
        constexpr bool operator==(const iterator& other) const noexcept {
            return ptr_ == other.ptr_;
        }
        
        constexpr bool operator!=(const iterator& other) const noexcept {
            return ptr_ != other.ptr_;
        }
    };

    constexpr small_unordered_set() noexcept = default;
    
    ~small_unordered_set() noexcept = default;
    
    constexpr iterator begin() noexcept {
        if (size_ == 0) {
            return end();
        }
        return iterator(&data_[0], *this, 0);
    }
    
    constexpr iterator end() noexcept {
        return iterator(nullptr, *this, size_);
    }
    
    constexpr bool empty() const noexcept {
        return size_ == 0;
    }
    
    constexpr size_t size() const noexcept {
        return size_;
    }
    
    constexpr size_t max_size() const noexcept {
        return MAX_SIZE;
    }
    
    constexpr void clear() noexcept {
        size_ = 0;
    }
    
    constexpr bool insert(const T& value) noexcept {
        const T* existing = find_element(value);
        if (existing) {
            return false;
        }
        
        if (size_ >= MAX_SIZE) {
            error(np_debug_info, std::format("small_unordered_set::insert: max_size reached, current size: {}", size_));
            return false;
        }
        
        data_[size_] = value;
        size_++;
        return true;
    }
    
    constexpr bool erase(const T& value) noexcept {
        const T* existing = find_element(value);
        if (!existing) {
            return false;
        }
        
        size_t index = existing - data_;
        
        for (size_t i = index; i < size_ - 1; ++i) {
            data_[i] = data_[i + 1];
        }
        
        size_--;
        return true;
    }
    
    constexpr iterator find(const T& value) noexcept {
        T* found = find_element(value);
        if (found) {
            return iterator(found, this, found - data_);
        }
        return end();
    }
    
    constexpr size_t count(const T& value) const noexcept {
        return find_element(value) ? 1 : 0;
    }
    
    constexpr bool contains(const T& value) const noexcept {
        return find_element(value) != nullptr;
    }
    
    constexpr bool operator==(const small_unordered_set& other) const noexcept {
        if (size_ != other.size_) {
            return false;
        }
        
        for (size_t i = 0; i < size_; ++i) {
            if (!other.contains(data_[i])) {
                return false;
            }
        }
        
        return true;
    }
    
    constexpr bool operator!=(const small_unordered_set& other) const noexcept {
        return !(*this == other);
    }
};

} // namespace grapology

#endif // GRAPOLOGY_SMALL_SET_HPP