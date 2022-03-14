#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.
  BBox bbox;

  int num_primitives = 0;
  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
    num_primitives++;
  }
  
  BVHNode *node = new BVHNode(bbox);

  std::vector<Primitive *>::iterator bound;
 
  if(num_primitives <= max_leaf_size) {
    node->l = NULL;
    node->r = NULL;
    node->start = start;
    node->end = end;

  } else {
    int ax = 0; 

    Vector3D extent = node->bb.extent;

    if (max(max(extent.x, extent.y), extent.z) == extent.y) {
      ax = 1;
    }
    if (max(max(extent.x, extent.y), extent.z) == extent.z) {
      ax = 2;
    }
    
    double spl = 0.0;
 
    for (auto p = start; p != end; p++) {
      spl += (*p)->get_bbox().centroid()[ax];
    }
    spl /= num_primitives;
    bound = std::partition(start, end, [spl, ax](Primitive * a) -> bool { return a->get_bbox().centroid()[ax] < spl; });

    //prevent infinite recursion
    if (start == bound) {
      bound++;
    }

    if (end == bound) {
      bound--;
    }

    node->l = construct_bvh(start, bound, max_leaf_size);
    node->r = construct_bvh(bound, end, max_leaf_size);
  return node;
  }

}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.



   double t0 = ray.min_t;
  double t1 = ray.max_t;
  bool didIntersect = node->bb.intersect(ray, t0, t1);

  if (t0 < ray.min_t || t1 > ray.max_t) {
    return false;
  }

  if (didIntersect) {
    if (node->isLeaf()) {
      for (auto p = node->start; p != node->end; node++) {
        auto deref = *p;
        if(deref->has_intersection(ray)) {
          return true;
        }
      }
      return false;
    } else {
      return has_intersection(ray, node->l) | has_intersection(ray, node->r);
    }
  } else {
    return false;
  }


}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.


  double t0 = ray.min_t;
  double t1 = ray.max_t;

  bool didIntersect = node->bb.intersect(ray, t0, t1);

  std::cout << didIntersect << endl;

  // if (t0 < ray.min_t || t1 > ray.max_t) {
  //   return false;
  // }

  bool didHitPrim = false;
  if (didIntersect) {
    if (node->isLeaf()) {
      for (auto p = node->start; p != node->end; p++) {
        total_isects++;
        didHitPrim = didHitPrim | (*p)->intersect(ray, i);
      }
      return didHitPrim;
    } else {
      bool left = intersect(ray, i, node->l);
      bool right = intersect(ray, i, node->r);
      return left | right;
    } 
  } else {
    return false;
  }


}

} // namespace SceneObjects
} // namespace CGL
