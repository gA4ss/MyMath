#ifndef MYMATH_AABOXKDTREE2D_HPP
#define MYMATH_AABOXKDTREE2D_HPP

#include <limits>
#include <memory>
#include <vector>
#include <cassert>
#include <algorithm>

#include "math_utils.hpp"
#include "aabox2d.hpp"

namespace mypilot {
namespace mymath {

// 轴对齐盒子KD树的参数
struct AABoxKDTreeParams {
  int max_depth = -1;                   // kdtree最大深度
  int max_leaf_size = -1;               // 在一个叶子节点的最大项目数量
  double max_leaf_dimension = -1.0;     // 叶子节点的最大维度
};

// 轴对齐盒子KD树节点
template <class ObjectType>
class AABoxKDTree2dNode {
public:
  using ObjectPtr = const ObjectType *;

  // kdtree节点的深度
  AABoxKDTree2dNode(const std::vector<ObjectPtr>& objects,
                    const AABoxKDTreeParams& params, int depth) : 
    _depth(depth) {
    assert(!objects.empty());

    compute_boundary(objects);
    compute_partition();

    if (split_to_subnodes(objects, params)) {
      std::vector<ObjectPtr> left_subnode_objects;
      std::vector<ObjectPtr> right_subnode_objects;
      partition_objects(objects, &left_subnode_objects, &right_subnode_objects);

      // 切分子节点
      if (!left_subnode_objects.empty()) {
        _left_subnode.reset(new AABoxKDTree2dNode<ObjectType>(
          left_subnode_objects, params, depth + 1));
      }
      if (!right_subnode_objects.empty()) {
        _right_subnode.reset(new AABoxKDTree2dNode<ObjectType>(
          right_subnode_objects, params, depth + 1));
      }
    } else {
      init_objects(objects);
    }
  }

  // 通过以节点'point'为根的KD树，获取最接近目标点的对象。
  ObjectPtr get_nearest_object(const Vec2d& point) const {
    ObjectPtr nearest_object = nullptr;
    double min_distance_sqr = std::numeric_limits<double>::infinity();
    get_nearest_object_internal(point, &min_distance_sqr, &nearest_object);
    return nearest_object;
  }

  // 搜索以'point'为中心，'distance'为距离的范围，获取到点距离之内的对象。
  std::vector<ObjectPtr> get_objects(const Vec2d& point,
                                     const double distance) const {
    std::vector<ObjectPtr> result_objects;
    get_objects_internal(point, distance, square(distance), &result_objects);
    return result_objects;
  }

  // 获取轴对齐盒子
  AABox2d get_bounding_box() const {
    return AABox2d({_min_x, _min_y}, {_max_x, _max_y});
  }

private:
  void init_objects(const std::vector<ObjectPtr>& objects) {
    _num_objects = objects.size();
    _objects_sorted_by_min = objects;
    objects_sorted_by_max_ = objects;
    std::sort(_objects_sorted_by_min.begin(), _objects_sorted_by_min.end(),
              [&](ObjectPtr obj1, ObjectPtr obj2) {
                return _partition == PARTITION_X
                           ? obj1->aabox().min_x() < obj2->aabox().min_x()
                           : obj1->aabox().min_y() < obj2->aabox().min_y();
              });
    std::sort(objects_sorted_by_max_.begin(), objects_sorted_by_max_.end(),
              [&](ObjectPtr obj1, ObjectPtr obj2) {
                return _partition == PARTITION_X
                           ? obj1->aabox().max_x() > obj2->aabox().max_x()
                           : obj1->aabox().max_y() > obj2->aabox().max_y();
              });
    _objects_sorted_by_min_bound.reserve(_num_objects);
    for (ObjectPtr object : _objects_sorted_by_min) {
      _objects_sorted_by_min_bound.push_back(_partition == PARTITION_X
                                                 ? object->aabox().min_x()
                                                 : object->aabox().min_y());
    }
    _objects_sorted_by_max_bound.reserve(_num_objects);
    for (ObjectPtr object : objects_sorted_by_max_) {
      _objects_sorted_by_max_bound.push_back(_partition == PARTITION_X
                                                 ? object->aabox().max_x()
                                                 : object->aabox().max_y());
    }
  }

  bool split_to_subnodes(const std::vector<ObjectPtr>& objects,
                         const AABoxKDTreeParams& params) {
    if (params.max_depth >= 0 && _depth >= params.max_depth) {
      return false;
    }
    if (static_cast<int>(objects.size()) <= std::max(1, params.max_leaf_size)) {
      return false;
    }
    if (params.max_leaf_dimension >= 0.0 &&
        std::max(_max_x - _min_x, _max_y - _min_y) <=
          params.max_leaf_dimension) {
      return false;
    }
    return true;
  }

  double lower_distance_square_to_point(const Vec2d& point) const {
    double dx = 0.0;
    if (point.x() < _min_x) {
      dx = _min_x - point.x();
    } else if (point.x() > _max_x) {
      dx = point.x() - _max_x;
    }
    double dy = 0.0;
    if (point.y() < _min_y) {
      dy = _min_y - point.y();
    } else if (point.y() > _max_y) {
      dy = point.y() - _max_y;
    }
    return dx * dx + dy * dy;
  }

  double upper_distance_square_to_point(const Vec2d& point) const {
    const double dx =
      (point.x() > _mid_x ? (point.x() - _min_x) : (point.x() - _max_x));
    const double dy =
      (point.y() > _mid_y ? (point.y() - _min_y) : (point.y() - _max_y));
    return dx * dx + dy * dy;
  }

  void get_all_objects(std::vector<ObjectPtr>* const result_objects) const {
    result_objects->insert(result_objects->end(),
                           _objects_sorted_by_min.begin(),
                           _objects_sorted_by_min.end());
    if (_left_subnode != nullptr) {
      _left_subnode->get_all_objects(result_objects);
    }
    if (_right_subnode != nullptr) {
      _right_subnode->get_all_objects(result_objects);
    }
  }

  void get_objects_internal(const Vec2d& point, const double distance,
                            const double distance_sqr,
                            std::vector<ObjectPtr>* const result_objects) const {
    if (lower_distance_square_to_point(point) > distance_sqr) {
      return;
    }
    if (upper_distance_square_to_point(point) <= distance_sqr) {
      get_all_objects(result_objects);
      return;
    }
    const double pvalue = (_partition == PARTITION_X ? point.x() : point.y());
    if (pvalue < _partition_position) {
      const double limit = pvalue + distance;
      for (int i = 0; i < _num_objects; ++i) {
        if (_objects_sorted_by_min_bound[i] > limit) {
          break;
        }
        ObjectPtr object = _objects_sorted_by_min[i];
        if (object->distance_square_to(point) <= distance_sqr) {
          result_objects->push_back(object);
        }
      }
    } else {
      const double limit = pvalue - distance;
      for (int i = 0; i < _num_objects; ++i) {
        if (_objects_sorted_by_max_bound[i] < limit) {
          break;
        }
        ObjectPtr object = objects_sorted_by_max_[i];
        if (object->distance_square_to(point) <= distance_sqr) {
          result_objects->push_back(object);
        }
      }
    }
    if (_left_subnode != nullptr) {
      _left_subnode->get_objects_internal(point, distance, distance_sqr,
                                          result_objects);
    }
    if (_right_subnode != nullptr) {
      _right_subnode->get_objects_internal(point, distance, distance_sqr,
                                           result_objects);
    }
  }

  void get_nearest_object_internal(const Vec2d& point,
                                   double* const min_distance_sqr,
                                   ObjectPtr* const nearest_object) const {
    if (lower_distance_square_to_point(point) >= *min_distance_sqr - math_epsilon) {
      return;
    }
    const double pvalue = (_partition == PARTITION_X ? point.x() : point.y());
    const bool search_left_first = (pvalue < _partition_position);
    if (search_left_first) {
      if (_left_subnode != nullptr) {
        _left_subnode->get_nearest_object_internal(point, min_distance_sqr,
                                                   nearest_object);
      }
    } else {
      if (_right_subnode != nullptr) {
        _right_subnode->get_nearest_object_internal(point, min_distance_sqr,
                                                    nearest_object);
      }
    }
    if (*min_distance_sqr <= math_epsilon) {
      return;
    }

    if (search_left_first) {
      for (int i = 0; i < _num_objects; ++i) {
        const double bound = _objects_sorted_by_min_bound[i];
        if (bound > pvalue && square(bound - pvalue) > *min_distance_sqr) {
          break;
        }
        ObjectPtr object = _objects_sorted_by_min[i];
        const double distance_sqr = object->distance_square_to(point);
        if (distance_sqr < *min_distance_sqr) {
          *min_distance_sqr = distance_sqr;
          *nearest_object = object;
        }
      }
    } else {
      for (int i = 0; i < _num_objects; ++i) {
        const double bound = _objects_sorted_by_max_bound[i];
        if (bound < pvalue && square(bound - pvalue) > *min_distance_sqr) {
          break;
        }
        ObjectPtr object = objects_sorted_by_max_[i];
        const double distance_sqr = object->distance_square_to(point);
        if (distance_sqr < *min_distance_sqr) {
          *min_distance_sqr = distance_sqr;
          *nearest_object = object;
        }
      }
    }
    if (*min_distance_sqr <= math_epsilon) {
      return;
    }
    if (search_left_first) {
      if (_right_subnode != nullptr) {
        _right_subnode->get_nearest_object_internal(point, min_distance_sqr,
                                                    nearest_object);
      }
    } else {
      if (_left_subnode != nullptr) {
        _left_subnode->get_nearest_object_internal(point, min_distance_sqr,
                                                   nearest_object);
      }
    }
  }

  void compute_boundary(const std::vector<ObjectPtr>& objects) {
    _min_x = std::numeric_limits<double>::infinity();
    _min_y = std::numeric_limits<double>::infinity();
    _max_x = -std::numeric_limits<double>::infinity();
    _max_y = -std::numeric_limits<double>::infinity();
    for (ObjectPtr object : objects) {
      _min_x = std::fmin(_min_x, object->aabox().min_x());
      _max_x = std::fmax(_max_x, object->aabox().max_x());
      _min_y = std::fmin(_min_y, object->aabox().min_y());
      _max_y = std::fmax(_max_y, object->aabox().max_y());
    }
    _mid_x = (_min_x + _max_x) / 2.0;
    _mid_y = (_min_y + _max_y) / 2.0;
    assert(!std::isinf(_max_x) && 
           !std::isinf(_max_y) && 
           !std::isinf(_min_x) &&
           !std::isinf(_min_y));
    // the provided object box size is infinity
  }

  void compute_partition() {
    if (_max_x - _min_x >= _max_y - _min_y) {
      _partition = PARTITION_X;
      _partition_position = (_min_x + _max_x) / 2.0;
    } else {
      _partition = PARTITION_Y;
      _partition_position = (_min_y + _max_y) / 2.0;
    }
  }

  void partition_objects(const std::vector<ObjectPtr>& objects,
                         std::vector<ObjectPtr>* const left_subnode_objects,
                         std::vector<ObjectPtr>* const right_subnode_objects) {
    left_subnode_objects->clear();
    right_subnode_objects->clear();
    std::vector<ObjectPtr> other_objects;
    if (_partition == PARTITION_X) {
      for (ObjectPtr object : objects) {
        if (object->aabox().max_x() <= _partition_position) {
          left_subnode_objects->push_back(object);
        } else if (object->aabox().min_x() >= _partition_position) {
          right_subnode_objects->push_back(object);
        } else {
          other_objects.push_back(object);
        }
      }
    } else {
      for (ObjectPtr object : objects) {
        if (object->aabox().max_y() <= _partition_position) {
          left_subnode_objects->push_back(object);
        } else if (object->aabox().min_y() >= _partition_position) {
          right_subnode_objects->push_back(object);
        } else {
          other_objects.push_back(object);
        }
      }
    }
    init_objects(other_objects);
  }

private:
  int _num_objects = 0;
  std::vector<ObjectPtr> _objects_sorted_by_min;
  std::vector<ObjectPtr> objects_sorted_by_max_;
  std::vector<double> _objects_sorted_by_min_bound;
  std::vector<double> _objects_sorted_by_max_bound;
  int _depth = 0;

  // 边界
  double _min_x = 0.0;
  double _max_x = 0.0;
  double _min_y = 0.0;
  double _max_y = 0.0;
  double _mid_x = 0.0;
  double _mid_y = 0.0;

  enum Partition {
    PARTITION_X = 1,
    PARTITION_Y = 2,
  };
  Partition _partition = PARTITION_X;
  double _partition_position = 0.0;

  std::unique_ptr<AABoxKDTree2dNode<ObjectType>> _left_subnode = nullptr;
  std::unique_ptr<AABoxKDTree2dNode<ObjectType>> _right_subnode = nullptr;
};

// 轴对齐kdtree
template <class ObjectType>
class AABoxKDTree2d {
public:
  using ObjectPtr = const ObjectType *;

  AABoxKDTree2d(const std::vector<ObjectType>& objects,
                const AABoxKDTreeParams& params) {
    if (!objects.empty()) {
      std::vector<ObjectPtr> object_ptrs;
      for (const auto& object : objects) {
        object_ptrs.push_back(&object);
      }
      _root.reset(new AABoxKDTree2dNode<ObjectType>(object_ptrs, params, 0));
    }
  }

  // 获取离点'point'最近的对象
  ObjectPtr get_nearest_object(const Vec2d& point) const {
    return _root == nullptr ? nullptr : _root->get_nearest_object(point);
  }

  // 获取所有在树上的对象
  std::vector<ObjectPtr> get_objects(const Vec2d& point,
                                     const double distance) const {
    if (_root == nullptr) {
      return {};
    }
    return _root->get_objects(point, distance);
  }

  // 获取轴对齐盒子包含树里所有的对象
  AABox2d get_bounding_box() const {
    return _root == nullptr ? AABox2d() : _root->get_bounding_box();
  }

private:
  std::unique_ptr<AABoxKDTree2dNode<ObjectType>> _root = nullptr;
};

}}

#endif
