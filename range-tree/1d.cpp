#include <algorithm>
#include <iostream>
#include <iterator>
#include <queue>
#include <random>
#include <set>
#include <vector>

using std::vector;

struct Node {
  using NodePtr = Node*;

public:
  Node() = delete;
  Node(int val, NodePtr l = nullptr, NodePtr r = nullptr)
      : value(val), left(l), right(r), subtree_count(0) {}
  ~Node() {
    if (left)
      delete left;
    if (right)
      delete right;
  }
  Node(std::vector<int>::iterator beg, std::vector<int>::iterator end)
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
  }

  bool is_leaf() const { return !left && !right; }
  NodePtr rebuild() const {
    auto data = return_subtree_sorted();
    return new Node(data.begin(), data.end());
  }

  vector<int> return_subtree_sorted() const {
    if (!left && !right)
      return {value};

    vector<int> res;
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

  NodePtr insert(int x) {
    // leaf node
    if (!left && !right) {
      int smaller = std::min(x, value);
      int bigger = std::max(x, value);
      value = smaller;
      left = new Node(smaller);
      right = new Node(bigger);
      return this;
    }

    // internal node
    subtree_count += 1;

    if (x <= value) {
      left = left->insert(x);
      // update value of internal node
      value = std::max(x, value);
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
  int value;
  NodePtr left, right;
  int subtree_count;
};

class Tree {
  using NodePtr = Node*;

public:
  Tree() : root_(nullptr) {}
  void insert(int val) {
    if (!root_) {
      root_ = new Node(val);
      return;
    }
    root_ = root_->insert(val);
  }

  void print() const {
    if (!root_) {
      std::cout << "Empty tree!\n";
      return;
    }
    root_->print(0);
  }

  vector<int> between(int a, int b) const {
    if (!root_)
      return {};

    std::vector<int> out;

    NodePtr common = root_;
    while (!common->is_leaf()) {
      bool lpath_dir = a <= common->value;
      bool rpath_dir = b <= common->value;

      if (lpath_dir != rpath_dir)
        break;

      common = lpath_dir ? common->left : common->right;
    }

    if (common->is_leaf()) {
      if (int x = common->value; x > a && x < b)
        return {x};
      return {};
    }

    NodePtr lpath = common->left, rpath = common->right;

    while (!lpath->is_leaf()) {
      bool lpath_dir = a <= lpath->value;
      // going to the left
      if (lpath_dir) {
        auto temp = lpath->right->return_subtree_sorted();
        out.insert(out.end(), temp.begin(), temp.end());
      }
      lpath = lpath_dir ? lpath->left : lpath->right;
    }

    while (!rpath->is_leaf()) {
      bool rpath_dir = b <= rpath->value;
      if (!rpath_dir) {
        auto temp = rpath->left->return_subtree_sorted();
        out.insert(out.end(), temp.begin(), temp.end());
      }
      rpath = rpath_dir ? rpath->left : rpath->right;
    }

    if (lpath->value > a && lpath->value < b)
      out.push_back(lpath->value);
    if (rpath->value > a && rpath->value < b)
      out.push_back(rpath->value);

    return out;
  }

private:
  NodePtr root_;
};

vector<int> between(int a, int b, const std::vector<int>& vals) {
  vector<int> out;
  for (int x : vals)
    if (x > a && x < b) {
      out.push_back(x);
    }
  return out;
}

vector<int> between_alt(int a, int b, const Tree& t) {
  auto out = t.between(a, b);
  sort(out.begin(), out.end());
  return out;
}

int main() {
  srand(2104);
  Tree t;

  std::set<int> S;

  for (int i = 0; i < 20; ++i) {
    int r = rand() % 100;
    if (S.count(r) != 0)
      continue;
    S.insert(r);
    // std::cout << "Inserting " << r << "\n";
    t.insert(r);
    // std::cout << "Inserted successfully.\n";
    //  t.print();
  }

  std::vector<int> vals(S.begin(), S.end());

  std::cout << "Vals are:";
  for (int x : vals) {
    std::cout << " " << x;
  }
  std::cout << "\n";

  Node* test = new Node(vals.begin(), vals.end());
  test->print(0);

  for (int i = 0; i < 100; ++i) {
    int k = rand() % vals.size();
    int l = rand() % vals.size();

    if (k == l)
      continue;

    if (k > l)
      std::swap(k, l);

    int a = vals[k], b = vals[l];

    auto res1 = between(a, b, vals);
    auto res2 = between_alt(a, b, t);

    /*std::cout << "Between " << a << " and " << b << " is:";
    for (auto x : out) {
      std::cout << " " << x;
    }*/

    if (res1 != res2)
      std::cout << "Fail\n";
    /*else
      std::cout << "Succ\n";*/
  }

  std::cout << "Test:\n";
  Tree f;
  for (int i = 0; i < 100; ++i)
    f.insert(i);
  f.print();

  return 0;
}
