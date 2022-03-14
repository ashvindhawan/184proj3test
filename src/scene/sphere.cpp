#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.

  float a = dot(r.d, r.d);
  float b = dot(2*(r.o-o), r.d);
  float c = dot(r.o - o, r.o - o) - r2;

  float radicand = b*b - 4 * a * c;

  if(radicand < 0) {
    return false;
  }

  float t_plus = (-b + sqrt(radicand))/(2*a); 
  float t_minus = (-b - sqrt(radicand))/(2*a); 

  if (t_plus >= 0 || t_minus >= 0) {
    t1 = min(t_plus, t_minus);
    t2 = max(t_plus, t_minus);
    // std::cout <<"test t1: " << t1 << endl;
    // std::cout <<"test t2: " << t2 << endl;
    return true;
  } else {
    return false;
  }
}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1;
  double t2;
  if(test(r, t1, t2)) { //test function does equation evaluation, checks that t >= 0, and deals with undefined square roots
    if(t1>=r.min_t && t1<=r.max_t) {
      return true;
    }
    if(t2>=r.min_t && t2<=r.max_t) {
      return true;
    }
  }
  return false;
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.

  double t1;
  double t2;
  if(test(r, t1, t2)) { //test function does equation evaluation, checks that t >= 0, and deals with undefined square roots
    if(t1>=r.min_t && t1<=r.max_t) {
      
      Vector3D int_point = r.o + t1 * r.d;
      Vector3D normal = (int_point-o).unit(); 
      r.max_t = t1;
      i->t = t1;
      i->n = normal;
      i->primitive = this;
      i->bsdf = get_bsdf();
      return true;

    } else if(t2>=r.min_t && t2<=r.max_t) {
      Vector3D int_point = r.o + t2 * r.d;
      Vector3D normal = (int_point-o).unit(); 
      r.max_t = t2;
      i->t = t2;
      i->n = normal;
      i->primitive = this;
      i->bsdf = get_bsdf();
      return true;
    }
  }
  return false;
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
