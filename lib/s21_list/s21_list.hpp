#pragma once

#include <initializer_list>
#include <iostream>
#include <limits>

/*
  synopsis:
  https://en.cppreference.com/w/cpp/container/list
  https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/n4835.pdf#page=824&zoom=100,85,710

  copy-and-swap:
  https://ru.wikipedia.org/wiki/Copy-and-swap

  node:
  https://www.codesdope.com/blog/article/c-linked-lists-in-c-singly-linked-list/
*/

namespace s21 {
template <typename T> class List {
  class ListNode;
  class ListIterator;
  class ListConstIterator;

public:
  // ========================== *List Member type* ==========================
  using value_type = T;
  using reference = value_type &;
  using const_reference = const value_type &;
  using iterator = ListIterator;
  using const_iterator = ListConstIterator;
  using size_type = std::size_t;

  // ========================== *List Functions* ==========================
  List() : head_(nullptr), tail_(nullptr), size_(0) {
    // создаем end_ заранее, сама нода показывает кол-во элементов в списке
    end_ = new ListNode(value_type{});
  }

  explicit List(size_type n) : head_(nullptr), tail_(nullptr), size_(0) {
    end_ = new ListNode(size_);
    while (n > 0) {
      push_back(value_type());
      n--;
    }
  }

  List(std::initializer_list<value_type> const &items)
      : head_(nullptr), tail_(nullptr), size_(0) {
    end_ = new ListNode(size_);
    for (auto it = items.begin(); it != items.end(); ++it) {
      push_back(*it); // const_iterator
    }
  }

  List(const List &other) : head_(nullptr), tail_(nullptr), size_(0) {
    end_ = new ListNode(size_);
    for (auto it = other.cbegin(); it != other.cend(); ++it) {
      push_back(*it); // const_iterator
    }
  }

  List(List &&other) noexcept
      : head_(nullptr), tail_(nullptr), end_(nullptr), size_(0) {
    this->Swap(other);
  }

  ~List() {
    Clear();
    if (end_) {
      delete end_;
    }
  }

  // copy
  List &operator=(const List &other) {
    if (this != &other) {
      List tmp(other); // create copy other
      this->Swap(tmp); // same *this = std::move(tmp);
    }
    return *this;
  }
  // move
  List &operator=(List &&other) noexcept {
    if (this != &other) {
      this->Swap(other);
    }
    return *this;
  }

  // ========================== *List Element access* ==========================
  reference Front() {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    return *begin();
  }
  const_reference Front() const {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    return *cbegin();
  }
  reference Back() {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    return tail_->value_;
  }
  const_reference Back() const {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    return tail_->value_;
  }

  // ========================== *List Iterators* ==========================
  iterator begin() noexcept {
    if (Empty()) {
      // *begin = 0  =>  return *iterator(end_) = 0
      return iterator(end_);
    } else {
      return iterator(head_);
    }
  }
  const_iterator begin() const noexcept {
    if (Empty()) {
      return const_iterator(end_);
    } else {
      return const_iterator(head_);
    }
  }
  const_iterator cbegin() const noexcept {
    if (Empty()) {
      return const_iterator(end_);
    } else {
      return const_iterator(head_);
    }
  }

  // return end_, even if list empty. same (begin() + size())
  iterator end() noexcept { return iterator(end_); }
  const_iterator end() const noexcept { return const_iterator(end_); }
  const_iterator cend() const noexcept { return const_iterator(end_); }

  // ========================== *List Capacity* ==========================
  bool Empty() const noexcept { return (size_ == 0); }
  size_type Size() const noexcept { return size_; }

  size_type MaxSize() const noexcept {
    return (std::numeric_limits<size_type>::max()) / sizeof(ListNode);
    // if linux: max() / 2
    // if mac: max()
  }

  void Clear() {
    while (!Empty()) {
      PopBack();
    }
  }

  iterator Insert(const_iterator pos, const_reference value) {
    if (Empty()) {
      PushBack(value);
      return iterator(tail_);
    }
    if (pos == const_iterator(head_)) {
      PushFront(value);
      return iterator(head_);
    }
    if (pos == const_iterator(end_)) {
      PushBack(value);
      return iterator(tail_);
    }
    // create ptr for left and right node
    ListNode *right = pos.ptr_;
    ListNode *left = pos.ptr_->prev_;
    // create new node for list
    ListNode *tmp = new ListNode(value);
    tmp->prev_ = left;
    tmp->next_ = right;
    // change other nodes
    size_++;
    AddEnd(); // end_->value_ = size_;
    left->next_ = tmp;
    right->prev_ = tmp;
    return iterator(tmp);
  }

  void Erase(iterator pos) {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    if (pos == iterator(end_)) {
      throw std::out_of_range("ptr from list, err in original");
    }
    if (pos == iterator(head_)) {
      PopFront();
    } else if (pos == iterator(tail_)) {
      PopBack();
    } else {
      // create ptr for current, left and right node
      ListNode *tmp = pos.ptr_;
      ListNode *right = pos.ptr_->next_;
      ListNode *left = pos.ptr_->prev_;
      left->next_ = right;
      right->prev_ = left;
      delete tmp;
      size_--;
      AddEnd();
    }
  }

  void PushBack(const_reference value) {
    ListNode *tmp = new ListNode(value);
    tmp->prev_ = tail_;
    if (Empty()) {
      head_ = tmp;
      tail_ = tmp;
    } else {
      tail_->next_ = tmp;
      tail_ = tmp;
    }
    size_++;
    AddEnd();
  }

  void PopBack() {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    if (size_ == 1) {
      delete head_;
      head_ = nullptr;
      tail_ = nullptr;
      size_ = 0;
    } else {
      ListNode *tmp = tail_->prev_; // предпоследний элемент в листе
      tmp->next_ = nullptr;
      delete tail_;
      tail_ = tmp;
      size_--;
    }
    AddEnd();
  }

  void PushFront(const_reference value) {
    ListNode *tmp = new ListNode(value);
    tmp->next_ = head_;
    if (Empty()) {
      head_ = tmp;
      tail_ = tmp;
    } else {
      head_->prev_ = tmp;
      head_ = tmp;
    }
    size_++;
    AddEnd();
  }

  void PopFront() {
    if (Empty()) {
      throw std::out_of_range("list is empty");
    }
    if (size_ == 1) {
      delete head_;
      head_ = nullptr;
      tail_ = nullptr;
      size_ = 0;
    } else {
      ListNode *tmp = head_->next_; // 2ой элемент в листе
      tmp->prev_ = nullptr;
      delete head_;
      head_ = tmp;
      size_--;
    }
    AddEnd();
  }

  void Swap(List &other) noexcept {
    std::swap(this->head_, other.head_);
    std::swap(this->tail_, other.tail_);
    std::swap(this->end_, other.end_);
    std::swap(this->size_, other.size_);
  }

  void Merge(List &other) {
    // если this пустой
    if (Empty()) {
      *this = std::move(other); // через оператор= перемещения, тк в оригинале
                                // other стирается
    } else if (!other.Empty()) {
      // если other пустой, то ничего не происходит, иначе мерджим
      const_iterator it_this = const_iterator(head_);
      const_iterator it_other = const_iterator(other.head_);
      // идем по листу this
      while (it_this != const_iterator(end_)) {
        if (*it_this > *it_other) {
          Insert(it_this, *it_other);
          ++it_other;
        } else {
          ++it_this;
        }
      }
    }
    other.Clear();
  }

  void Splice(const_iterator pos, List &other) {
    if (Empty()) {
      *this = std::move(other);
    } else if (!other.Empty()) {
      const_iterator it = const_iterator(other.head_);
      while (it != const_iterator(other.end_)) {
        Insert(pos, *it);
        ++it;
      }
      other.Clear();
    }
  }

  void Reverse() {
    if (!Empty() && (size_ != 1)) {
      iterator left = iterator(head_);
      iterator right = iterator(tail_);
      size_type i = 0;
      while (i < size_ / 2) {
        std::swap(*left, *right);
        ++left;
        --right;
        i++;
      }
    }
  }

  void Unique() {
    if (!Empty() && (size_ != 1)) {
      size_type i = 0;
      iterator left = iterator(head_);
      iterator right = iterator(head_->next_);
      while (i < size_) {
        if (*left == *right) {
          Erase(right);
          right = left;
          ++right;
        } else {
          ++left;
          ++right;
        }
        i++;
      }
    }
  }

  void Sort() {
    if (!Empty() && (size_ != 1)) {
      iterator it_j = iterator(head_);
      iterator it_j_next = iterator(head_->next_);
      for (size_type i = 0; i < Size() - 1; i++) {
        // return pointers to the beginning
        it_j = iterator(head_);
        it_j_next = iterator(head_->next_);
        for (size_type j = 0; j < (Size() - i - 1); j++) {
          if (*it_j > *it_j_next) {
            std::swap(*it_j, *it_j_next);
          }
          ++it_j;
          ++it_j_next;
        }
      }
    }
  }
  // ========================== *Bonus* ==========================
  template <class... Args>
  iterator InsertMany(const_iterator pos, Args &&...args) {
    iterator it = nullptr;
    for (const auto &arg : {args...}) {
      it = Insert(pos, arg);
    }
    return it;
  }

  template <class... Args> void InsertManyBack(Args &&...args) {
    Insert_many(const_iterator(end_), args...);
  }

  template <class... Args> void InsertManyFront(Args &&...args) {
    Insert_many(const_iterator(head_), args...);
  }

  // ========================== *Help functions* ==========================
  void s21PrintfList() {
    std::cout << "begin:" << *begin() << '\n';
    std::cout << "end:" << *end() << "\n";
    std::cout << "my_list: ";
    for (iterator it = begin(); it != end(); ++it) {
      std::cout << *it << " ";
    }
    std::cout << "\n";
  }

  void AddEnd() {
    end_->prev_ = tail_;
    end_->next_ = head_;
    end_->value_ = value_type{};
    // если у нас есть хоть 1 элемент связываем с end_
    // проверка на пустоту иначе сега
    if (!Empty()) {
      tail_->next_ = end_;
      head_->prev_ = end_;
    }
  }

private:
  // ========================== *Node* ==========================
  class ListNode {
  public:
    value_type value_;
    ListNode *next_;
    ListNode *prev_;

    ListNode(const_reference value) {
      this->value_ = value;
      this->next_ = nullptr;
      this->prev_ = nullptr;
    }
  };

  // ========================== *Iteraror* ==========================
  class ListIterator {
  public:
    ListNode *ptr_ = nullptr;
    ListIterator() { ptr_ = nullptr; }
    ListIterator(ListNode *ptr) { ptr_ = ptr; } // update ptr
    ~ListIterator() { ptr_ = nullptr; }

    reference operator*() {
      if (ptr_) {
        return ptr_->value_;
      }
      throw std::out_of_range("error");
    }

    bool operator==(const iterator &other) { return (ptr_ == other.ptr_); }

    bool operator!=(const iterator &other) { return (ptr_ != other.ptr_); }

    iterator &operator++() {
      ptr_ = ptr_->next_;
      return *this;
    }
    iterator &operator--() {
      ptr_ = ptr_->prev_;
      return *this;
    }
  };

  // ========================== *Const Iteraror* ==========================
  class ListConstIterator {
  public:
    ListNode *ptr_;
    ListConstIterator() { ptr_ = nullptr; }
    ListConstIterator(ListNode *ptr) { ptr_ = ptr; }
    ~ListConstIterator() { ptr_ = nullptr; }

    const_reference operator*() {
      if (ptr_) {
        return ptr_->value_;
      }
      throw std::out_of_range("error");
    }

    bool operator==(const const_iterator &other) {
      return (ptr_ == other.ptr_);
    }

    bool operator!=(const const_iterator &other) {
      return (ptr_ != other.ptr_);
    }

    const_iterator &operator++() {
      ptr_ = ptr_->next_;
      return *this;
    }

    const_iterator &operator--() {
      ptr_ = ptr_->prev_;
      return *this;
    }
  };

  // ========================== *List* ==========================
  ListNode *head_;
  ListNode *tail_;
  ListNode *end_;  // доп узел в конце листа
  size_type size_; // кол-во элементов в листе

}; // list
} // namespace s21
