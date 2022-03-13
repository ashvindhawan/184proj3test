#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

Vector3D Triangle::mullerTrumbore(const Ray &r) const {
  Vector3D e1 = p2 - p1;
  Vector3D e2 = p3 - p1;
  Vector3D s = r.o - p1;
  Vector3D s1 = cross(r.d, e2);
  Vector3D s2 = cross(s, e1);

  Vector3D answer(dot(s2,e2), dot(s1,s), dot(s2,r.d)); 
  answer /= (dot(s1,e1));
  return answer;
}

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.

  Vector3D answer = mullerTrumbore(r);

  if(answer[0] < r.min_t || answer[0] > r.max_t) {
    return false;
  }
  if (answer[1] < 0 || answer[1] > 1) {
    return false;
  }
  if (answer[2] < 0 || answer[2] > 1) {
    return false;
  }
  if (1 - answer[1] - answer[2] < 0 || 1 - answer[1] - answer[2] > 1) {
    return false;
  }

  // at this point we know it's true
  r.max_t = answer[0];
  return true;

}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly

    Vector3D answer = mullerTrumbore(r);

  if(answer[0] < r.min_t || answer[0] > r.max_t) {
    return false;
  }
  if (answer[1] < 0 || answer[1] > 1) {
    return false;
  }
  if (answer[2] < 0 || answer[2] > 1) {
    return false;
  }
  if (1 - answer[1] - answer[2] < 0) {
    return false;
  }

// at this point we know it's true
  r.max_t = answer[0];
  isect->t = answer[0];
  isect->n = answer[1] * n2 + answer[2] * n3 + (1-answer[1]-answer[2]) * n1;
  isect->primitive = this;
  isect->bsdf = get_bsdf();
  return true;
}



void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
