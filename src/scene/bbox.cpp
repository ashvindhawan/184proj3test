#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
    double entry_x = (min.x - r.o.x) / r.d.x;
    double entry_y = (min.y - r.o.y) / r.d.y;
    double entry_z = (min.z - r.o.z) / r.d.z;
    double exit_x = (max.x - r.o.x) / r.d.x;
    double exit_y = (max.y - r.o.y) / r.d.y;
    double exit_z = (max.z - r.o.z) / r.d.z;

    if (entry_x > exit_x) {
      std::swap(entry_x, exit_x);
    }
    if (entry_y > exit_y) {
      std::swap(entry_y, exit_y);
    }
    if (entry_z > exit_z) {
      std::swap(entry_z, exit_z);
    }

    double new_t0 = std::max(entry_x, std::max(entry_y, entry_z));
    double new_t1 = std::min(exit_x, std::min(exit_y, exit_z));

    // if (new_t0 < t0 || new_t1 > t1 || new_t1 < new_t0) {
    //   return false;
    // }

    if(new_t1 < new_t0) {
      return false;
    }
    // std::cout << "old_t0: " << t0 << "  old_t1: " <<  t1 << std::endl;
    // std::cout << "new_t0: " << new_t0 << "  new_t1: " <<  new_t1 << std::endl;

    t0 = new_t0;
    t1 = new_t1;
  return true;

}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
