#pragma once

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using std::string;
using std::vector;

/*
class Point {
public:
  Point() = default;
  Point(double x, double y, double z) : x_(x), y_(y), z_(z) {}

  template <std::size_t Index> std::tuple_element_t<Index, Point> get() const {
    if constexpr (Index == 0)
      return x_;
    if constexpr (Index == 1)
      return y_;
    if constexpr (Index == 2)
      return z_;
  }

  friend std::ostream& operator<<(std::ostream& os, const Point& x) {
    os << '[' << x.x_ << ',' << x.y_ << ',' << x.z_ << ']';
    return os;
  }

  friend bool operator==(const Point& a, const Point& b) {
    return (a.x_ == b.x_) && (a.y_ == b.y_) && (a.z_ == b.z_);
  }

private:
  double x_;
  double y_;
  double z_;
};

namespace std {
template <> struct tuple_size<::Point> { static constexpr size_t value = 3; };

template <> struct tuple_element<0, ::Point> { using type = double; };

template <> struct tuple_element<1, ::Point> { using type = double; };

template <> struct tuple_element<2, ::Point> { using type = double; };

} // namespace std
*/

template <typename DataType, size_t dim> class Tree;
template <typename DataType, size_t dim> struct Node;

// node in the last dimension uses this as TreeType
class TreeEmpty {};

// holds value, left and right child but also poiner to higher dimension tree
template <typename DataType, size_t dim> struct Node {
  using NodeType = Node<DataType, dim>;
  using NodePtr = NodeType*;
  static constexpr bool is_last_dim =
      (dim + 1) == std::tuple_size<DataType>::value;
  using TreeType = typename std::conditional<is_last_dim, TreeEmpty,
                                             Tree<DataType, dim + 1>>::type;
  using IndexType = typename std::tuple_element<dim, DataType>::type;

  static constexpr auto dim_getter = [](const DataType& x) -> IndexType {
    return x.template get<dim>();
  };

public:
  Node() = delete;
  Node(DataType val, NodePtr l = nullptr, NodePtr r = nullptr)
      : value(val), left(l), right(r), tree(), subtree_count(0) {
    if constexpr (!is_last_dim)
      tree.insert(val);
  }
  Node(Node&& o) noexcept {
    left = right = nullptr;
    subtree_count = 0;
    std::swap(value, o.value);
    std::swap(subtree_count, o.subtree_count);
    std::swap(tree, o.tree);
    std::swap(left, o.left);
    std::swap(right, o.right);
  }
  ~Node() {
    if (left)
      delete left;
    if (right)
      delete right;
    left = right = nullptr;
  }
  Node(typename std::vector<DataType>::iterator beg,
       typename std::vector<DataType>::iterator end)
      : left(nullptr), right(nullptr), tree() {

    auto count = std::distance(beg, end);
    if (count == 1) {
      value = *beg;
      subtree_count = 1;
      if constexpr (!is_last_dim) {
        tree = TreeType(beg, end);
      }
      return;
    }
    auto mid = beg + count / 2;

    value = *(mid - 1);
    left = new NodeType(beg, mid);
    right = new NodeType(mid, end);
    subtree_count = count;
    if constexpr (!is_last_dim) {
      auto cmp = [](const auto &a, const auto &b) {
        return a.template get<dim + 1>() < b.template get<dim + 1>();
      };
      std::inplace_merge(beg, mid, end, cmp);
      std::vector<DataType> out(beg, end);
      tree = TreeType(out.begin(), out.end());
    }
  }

  bool is_leaf() const { return !left && !right; }
  NodePtr rebuild() const {
    std::vector<DataType> data;
    data.reserve(subtree_count);
    return_subtree_sorted(data);
    return new NodeType(data.begin(), data.end());
  }

  void return_subtree_sorted(std::vector<DataType>& out) const {
    if (is_leaf())
      out.push_back(value);

    if (left)
      left->return_subtree_sorted(out);
    if (right)
      right->return_subtree_sorted(out);
  }

  NodePtr insert(DataType x) {
    if (is_leaf()) {
      DataType smaller = (dim_getter(x) < dim_getter(value)) ? x : value;
      DataType bigger = (dim_getter(x) >= dim_getter(value)) ? x : value;
      value = smaller;
      left = new NodeType(smaller);
      right = new NodeType(bigger);

      if constexpr (!is_last_dim)
        tree.insert(x);

      return this;
    }

    // internal node
    subtree_count += 1;

    if constexpr (!is_last_dim)
      tree.insert(x);

    if (dim_getter(x) <= dim_getter(value)) {
      left = left->insert(x);
      // update value of internal node
      value = (dim_getter(x) >= dim_getter(value)) ? x : value;
    } else {
      right = right->insert(x);
    }

    auto threshold = 2 * subtree_count / 3;
    if (left->subtree_count > threshold || right->subtree_count > threshold) {
      auto new_node = rebuild();
      delete this;
      return new_node;
    }
    return this;
  }

public:
  DataType value;
  NodePtr left, right;
  TreeType tree;
  int subtree_count;
};

template <typename DataType, size_t dim = 0> class Tree {
  static constexpr size_t tuple_entries = std::tuple_size<DataType>::value;
  using NodeType = Node<DataType, dim>;
  using NodePtr = NodeType*;
  static constexpr bool is_last_dim =
      (dim + 1) == std::tuple_size<DataType>::value;

public:
  Tree() : root_(nullptr) {}
  Tree(typename std::vector<DataType>::iterator beg,
       typename std::vector<DataType>::iterator end)
      : root_(new NodeType(beg, end)) {}
  Tree(Tree&& o) noexcept {
    root_ = nullptr;
    std::swap(root_, o.root_);
  }
  ~Tree() {
    if (root_) {
      delete root_;
    }
    root_ = nullptr;
  }
  Tree& operator=(Tree&& other) {
    if (this != &other) {
      if (root_) {
        delete root_;
        root_ = nullptr;
      }
      std::swap(root_, other.root_);
    }
    return *this;
  }

  void insert(DataType val) {
    if (!root_) {
      root_ = new NodeType(val);
      return;
    }
    root_ = root_->insert(val);
  }

  vector<DataType> between(DataType a, DataType b) const {
    static constexpr auto dim_getter = [](const DataType& x) {
      return x.template get<dim>();
    };

    if (!root_)
      return {};

    NodePtr common = root_;
    while (!common->is_leaf()) {
      bool lpath_dir = dim_getter(a) <= dim_getter(common->value);
      bool rpath_dir = dim_getter(b) <= dim_getter(common->value);

      if (lpath_dir != rpath_dir)
        break;

      common = lpath_dir ? common->left : common->right;
    }

    std::vector<DataType> out;

    if (common->is_leaf()) {
      if (auto x = dim_getter(common->value);
          x > dim_getter(a) && x < dim_getter(b)) {
        if constexpr (is_last_dim)
          out.push_back(common->value);
        else
          out = common->tree.between(a, b);
      }
      return out;
    }

    NodePtr lpath = common->left, rpath = common->right;

    while (!lpath->is_leaf()) {
      bool lpath_dir = dim_getter(a) <= dim_getter(lpath->value);
      if (lpath_dir) {
        if constexpr (is_last_dim) {
          lpath->right->return_subtree_sorted(out);
        } else {
          std::vector<DataType> temp;
          temp = lpath->right->tree.between(a, b);
          out.insert(out.end(), temp.begin(), temp.end());
        }
      }
      lpath = lpath_dir ? lpath->left : lpath->right;
    }

    while (!rpath->is_leaf()) {
      bool rpath_dir = dim_getter(b) <= dim_getter(rpath->value);
      if (!rpath_dir) {
        if constexpr (is_last_dim) {
          rpath->left->return_subtree_sorted(out);
        } else {
          auto temp = rpath->left->tree.between(a, b);
          out.insert(out.end(), temp.begin(), temp.end());
        }
      }
      rpath = rpath_dir ? rpath->left : rpath->right;
    }

    if (dim_getter(lpath->value) > dim_getter(a) &&
        dim_getter(lpath->value) < dim_getter(b)) {
      if constexpr (is_last_dim)
        out.push_back(lpath->value);
      else {
        auto temp = lpath->tree.between(a, b);
        out.insert(out.end(), temp.begin(), temp.end());
      }
    }

    if (dim_getter(rpath->value) > dim_getter(a) &&
        dim_getter(rpath->value) < dim_getter(b)) {
      if constexpr (is_last_dim)
        out.push_back(rpath->value);
      else {
        auto temp = rpath->tree.between(a, b);
        out.insert(out.end(), temp.begin(), temp.end());
      }
    }

    return out;
  }

private:
  NodePtr root_;
};
/*
std::vector<Point> between(const std::vector<Point>& control, Point a,
                           Point b) {
  std::vector<Point> out;
  for (auto x : control) {
    if (x.get<0>() > a.get<0>() && x.get<0>() < b.get<0>() &&
        x.get<1>() > a.get<1>() && x.get<1>() < b.get<1>() &&
        x.get<2>() > a.get<2>() && x.get<2>() < b.get<2>())
      out.push_back(x);
  }
  return out;
}

double get_rand_double() {
  static std::default_random_engine rnd{2104};
  static std::uniform_real_distribution<double> dist;
  return dist(rnd);
}

void single_test(int test_points) {
  Tree<Point> t;
  std::vector<Point> control;

  for (int i = 0; i < test_points; ++i) {
    auto p = Point(get_rand_double(), get_rand_double(), get_rand_double());
    t.insert(p);
    control.push_back(p);
  }

  double x = get_rand_double(), y = get_rand_double(), z = get_rand_double();
  double del = get_rand_double();

  auto a = Point(x, y, z);
  auto b = Point(x + del, y + del, z + del);

  auto out = t.between(a, b);
  auto real = between(control, a, b);

  auto cmp = [](Point a, Point b) {
    if (a.get<0>() != b.get<0>())
      return a.get<0>() < b.get<0>();
    if (a.get<1>() != b.get<1>())
      return a.get<1>() < b.get<1>();
    return a.get<2>() < b.get<2>();
  };

  sort(out.begin(), out.end(), cmp);
  sort(real.begin(), real.end(), cmp);

  std::cout << out.size() << " vs " << real.size() << "\n";

  if (out != real)
    std::cout << "Fail!\n";
}

int main() {
  srand(2104);

  constexpr int kTestCases = 10'000;
  for (int i = 0; i < 5; ++i) {
    single_test(kTestCases);
  }

  constexpr size_t kBenchCases = 10'000;

  std::vector<Point> inserts;
  inserts.reserve(kBenchCases);
  for (size_t i = 0; i < kBenchCases; ++i) {
    inserts.push_back(
        Point(get_rand_double(), get_rand_double(), get_rand_double()));
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  Tree<Point> t;
  for (const Point& x : inserts)
    t.insert(x);
  auto t2 = std::chrono::high_resolution_clock::now();

  auto duration_insert =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  std::cout << "Duration of insert is " << duration_insert << "\n";

  auto cmp = [](Point a, Point b) {
    if (a.get<0>() != b.get<0>())
      return a.get<0>() < b.get<0>();
    if (a.get<1>() != b.get<1>())
      return a.get<1>() < b.get<1>();
    return a.get<2>() < b.get<2>();
  };

  sort(inserts.begin(), inserts.end(), cmp);

  std::vector<std::pair<Point, Point>> queries;
  for (size_t i = 0; i < kBenchCases; ++i) {
    auto a = get_rand_double(), b = get_rand_double(), c = get_rand_double();
    auto d = a + get_rand_double(), e = b + get_rand_double(),
         f = c + get_rand_double();
    queries.push_back({Point(a, b, c), Point(d, e, f)});
  }

  auto t3 = std::chrono::high_resolution_clock::now();
  size_t checksum = 0;
  for (const auto& [from, to] : queries) {
    auto out = t.between(from, to);
    checksum += out.size();
  }
  auto t4 = std::chrono::high_resolution_clock::now();

  std::cout << "Checksum is " << checksum << "\n";

  auto duration_find =
      std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
  std::cout << "Duration of find is " << duration_find << "\n";

  // std::cout << "Out: ";
  // for (auto x : out)
  //   std::cout << " " << x;
  // std::cout << "\n";

  // std::cout << "Real: ";
  // for (auto x : real)
  //   std::cout << " " << x;
  // std::cout << "\n";

  return 0;
}
*/

/*std::string rand_str() {
  std::string out;
  for (int i = 0; i <= rand() % 4; ++i) {
    out.push_back('a' + rand() % 20);
  }
  return out;
}*/
