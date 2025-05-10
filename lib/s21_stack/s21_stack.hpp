#pragma once

#include <iostream>
#include <s21_list.hpp>

namespace s21 {
template <typename T> class Stack {
public:
  // *Member type
  using value_type = T;
  using reference = T &;
  using const_reference = const T &;
  using size_type = size_t;

  // *Stack Member functions*
  Stack() : list_() {}

  explicit Stack(std::initializer_list<value_type> const &items) : list_() {
    for (const auto &item : items) {
      Push(item);
    }
  }

  Stack(const Stack &s) : list_(s.list_) {}

  Stack(Stack &&s) noexcept : list_(std::move(s.list_)) {}

  ~Stack() { list_.Clear(); }

  Stack &operator=(const Stack &s) {
    list_ = s.list_;
    return *this;
  }

  Stack &operator=(Stack &&s) noexcept {
    list_ = std::move(s.list_);
    return *this;
  }

  // *Stack Element access*
  const_reference Top() const { return list_.Front(); }

  // *Stack Capacity*
  bool Empty() const noexcept { return (Size() == 0); }
  size_type Size() const noexcept { return (list_.Size()); }

  // *Stack Modifiers*
  void Push(const_reference value) { list_.PushFront(value); }

  void Pop() { list_.PopFront(); }

  void Swap(Stack &other) noexcept { list_.Swap(other.list_); }

  // *Bonus
  template <class... Args> void InsertManyFront(Args &&...args) {
    for (const auto &arg : {args...}) {
      Push(arg);
    }
  }

  // *Help function
  void s21PrintfStack() {
    std::cout << "Stack: ";
    while (Size() > 0) {
      std::cout << "size:" << Size() << "value:" << Top() << "   ";
      Pop();
    }
    std::cout << "\n";
  }

private:
  List<value_type> list_;
};
} // namespace s21
