#pragma once
#include <iostream>

#include <s21_list.hpp>

namespace s21 {
template <typename T>

class Queue {
public:
  // *Member type
  using value_type = T;
  using reference = T &;
  using const_reference = const T &;
  using size_type = size_t;

  // *Queue Member functions*
  Queue() : list_() {}

  explicit Queue(std::initializer_list<value_type> const &items)
      : list_(items) {}

  Queue(const Queue &s) : list_(s.list_) {}

  Queue(Queue &&s) noexcept : list_(std::move(s.list_)) {}

  ~Queue() { list_.Clear(); }

  Queue &operator=(const Queue &s) {
    list_ = s.list_;
    return *this;
  }

  Queue &operator=(Queue &&s) noexcept {
    list_ = std::move(s.list_);
    return *this;
  }

  // *Queue Element access*
  const_reference Front() const { return list_.Front(); }
  const_reference Back() const { return list_.Back(); }

  // *Queue Capacity*
  bool Empty() const noexcept { return (Size() == 0); }
  size_type Size() const noexcept { return (list_.Size()); }
  // *Queue Modifiers*
  void Push(const_reference value) { list_.PushBack(value); }

  void Pop() { list_.PopFront(); }

  void Swap(Queue &other) noexcept { list_.Swap(other.list_); }

  // *Bonus
  template <class... Args> void InsertManyBack(Args &&...args) {
    list_.InsertManyBack(std::forward<Args>(args)...);
  }

  // *Help function
  void s21PrintfQueue() {
    std::cout << "Queue: ";
    while (Size() > 0) {
      std::cout << "size:" << Size() << "value:" << Front() << "   ";
      Pop();
    }
    std::cout << "\n";
  }

private:
  List<value_type> list_;
};
} // namespace s21
