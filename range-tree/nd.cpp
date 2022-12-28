#include <algorithm>
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

// std::pair<int, string> p;
// auto [x, y] = p;

class Point {
public:
  Point() = default;
  Point(double x, double y, double z) : x_(x), y_(y), z_(z) {}

  template <std::size_t Index> std::tuple_element_t<Index, Point>& get() {
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

template <typename DataType, size_t dim> class Tree;
template <typename DataType, size_t dim> struct Node;

class TreeEmpty {};

template <typename DataType, size_t dim> struct Node {
  using NodeType = Node<DataType, dim>;
  using NodePtr = NodeType*;
  static constexpr bool is_last_dim =
      (dim + 1) == std::tuple_size<DataType>::value;
  using TreeType = typename std::conditional<is_last_dim, TreeEmpty,
                                             Tree<DataType, dim + 1>>::type;
  using IndexType = typename std::tuple_element<dim, DataType>::type;

public:
  Node() = delete;
  Node(DataType val, NodePtr l = nullptr, NodePtr r = nullptr)
      : value(val), left(l), right(r), subtree_count(0), tree() {
    if constexpr (!is_last_dim)
      tree.insert(val);
  }
  ~Node() {
    if (left)
      delete left;
    if (right)
      delete right;
  }
  /*Node(typename std::vector<DataType>::iterator beg,
       typename std::vector<DataType>::iterator end)
      : left(nullptr), right(nullptr) {
    auto count = std::distance(beg, end);
    if (count == 1) {
      value = *beg;
      subtree_count = 1;
      return;
    }
    auto mid = beg + count / 2;

    value = *(mid - 1);
    left = new Node(beg, mid);
    right = new Node(mid, end);
    subtree_count = count;
  }*/

  bool is_leaf() const { return !left && !right; }
  /*NodePtr rebuild() const {
    auto data = return_subtree_sorted();
    return new Node(data.begin(), data.end());
  }*/

  vector<DataType> return_subtree_sorted() const {
    if (!left && !right)
      return {value};

    vector<DataType> res;
    if (left)
      res = left->return_subtree_sorted();
    if (right) {
      auto temp = right->return_subtree_sorted();
      res.insert(res.end(), temp.begin(), temp.end());
    }
    return res;
  }

  void print(int depth) const {
    for (int i = 0; i < depth; ++i)
      std::cout << "\t";

    // leaf
    if (!left && !right)
      std::cout << "L";

    std::cout << value << "\n";

    if (left)
      left->print(depth + 1);
    if (right)
      right->print(depth + 1);
  }

  // z P dostal x-tu suradnicu
  // P.template get<x>();

  NodePtr insert(DataType x) {
    // leaf node
    if (!left && !right) {
      DataType smaller =
          (x.template get<dim>() < value.template get<dim>()) ? x : value;
      DataType bigger =
          (x.template get<dim>() >= value.template get<dim>()) ? x : value;
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

    if (x.template get<dim>() <= value.template get<dim>()) {
      left = left->insert(x);
      // update value of internal node
      value = (x.template get<dim>() >= value.template get<dim>()) ? x : value;
    } else {
      right = right->insert(x);
    }

    /*auto threshold = 2 * subtree_count / 3;
    if (left->subtree_count > threshold || right->subtree_count > threshold) {
      auto new_node = rebuild();
      delete this;
      return new_node;
    }*/
    return this;
  }

public:
  DataType value;
  TreeType tree;
  NodePtr left, right;
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
  ~Tree() { delete root_; }

  void insert(DataType val) {
    if (!root_) {
      root_ = new NodeType(val);
      return;
    }
    root_ = root_->insert(val);
  }

  void print() const {
    if (!root_) {
      std::cout << "Empty tree!\n";
      return;
    }
    root_->print(/*dim=*/0);
  }

  vector<DataType> between(DataType a, DataType b) const {
    if (!root_)
      return {};

    NodePtr common = root_;
    while (!common->is_leaf()) {
      bool lpath_dir =
          a.template get<dim>() <= common->value.template get<dim>();
      bool rpath_dir =
          b.template get<dim>() <= common->value.template get<dim>();

      if (lpath_dir != rpath_dir)
        break;

      common = lpath_dir ? common->left : common->right;
    }

    std::vector<DataType> out;

    if (common->is_leaf()) {
      if (auto x = common->value.template get<dim>();
          x > a.template get<dim>() && x < b.template get<dim>()) {
        if constexpr (is_last_dim)
          out.push_back(common->value);
        else
          out = common->tree.between(a, b);
      }
      return out;
    }

    NodePtr lpath = common->left, rpath = common->right;

    while (!lpath->is_leaf()) {
      bool lpath_dir =
          a.template get<dim>() <= lpath->value.template get<dim>();
      if (lpath_dir) {
        std::vector<DataType> temp;
        if constexpr (is_last_dim) {
          temp = lpath->right->return_subtree_sorted();
        } else {
          temp = lpath->right->tree.between(a, b);
        }
        out.insert(out.end(), temp.begin(), temp.end());
      }
      lpath = lpath_dir ? lpath->left : lpath->right;
    }

    while (!rpath->is_leaf()) {
      bool rpath_dir =
          b.template get<dim>() <= rpath->value.template get<dim>();
      if (!rpath_dir) {
        std::vector<DataType> temp;
        if constexpr (is_last_dim) {
          temp = rpath->left->return_subtree_sorted();
        } else {
          temp = rpath->left->tree.between(a, b);
        }
        out.insert(out.end(), temp.begin(), temp.end());
      }
      rpath = rpath_dir ? rpath->left : rpath->right;
    }

    if (lpath->value.template get<dim>() > a.template get<dim>() &&
        lpath->value.template get<dim>() < b.template get<dim>()) {
      if constexpr (is_last_dim)
        out.push_back(lpath->value);
      else {
        auto temp = lpath->tree.between(a, b);
        out.insert(out.end(), temp.begin(), temp.end());
      }
    }

    if (rpath->value.template get<dim>() > a.template get<dim>() &&
        rpath->value.template get<dim>() < b.template get<dim>()) {
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

std::string rand_str() {
  std::string out;
  for (int i = 0; i <= rand() % 4; ++i) {
    out.push_back('a' + rand() % 20);
  }
  return out;
}

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

void single_test(int test_points) {

  static std::default_random_engine rnd{2104};
  static std::uniform_real_distribution<double> dist;

  Tree<Point> t;
  std::vector<Point> control;

  for (int i = 0; i < test_points; ++i) {
    auto p = Point(dist(rnd), dist(rnd), dist(rnd));
    t.insert(p);
    control.push_back(p);
    // std::cout << "Inserting " << p << "\n";
  }

  // t.print();

  double x = dist(rnd), y = dist(rnd), z = dist(rnd);
  double del = dist(rnd);

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

  std::cout << real.size() << "\n";

  if (out != real)
    std::cout << "Fail!\n";
}

int main() {
  srand(2104);

  constexpr int kTestCases = 10000;
  for (int i = 0; i < 100; ++i) {
    single_test(kTestCases);
  }

  /*std::cout << "Out: ";
  for (auto x : out)
    std::cout << " " << x;
  std::cout << "\n";

  std::cout << "Real: ";
  for (auto x : real)
    std::cout << " " << x;
  std::cout << "\n";*/

  return 0;
}
