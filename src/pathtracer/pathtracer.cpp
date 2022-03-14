#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // TODO (Part 3):
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.

    /* 
  For loop num_samples times:
    randomly sample a vector wj (starting at hit_p)
    call BVH accel intersect 
    Set Li = itntersect.bdsf.getemission()
    Computer Li dot BSDF::f(wo, wi)
    Divide by pdf
    sum and average results
  */

  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;
  Vector3D sample_vec; double pdf;
    for (int i = 0; i < num_samples; i++) {
      Vector3D Fr = isect.bsdf->sample_f(w_out, &sample_vec, &pdf); // this function gives us a random vector and calculates Fr
      sample_vec = o2w * sample_vec;
      Ray randomR = Ray(hit_p, sample_vec); // wasn't working unless we use the ray constructor??
      randomR.min_t = EPS_F;
      Intersection randomIsect;
      bool didIntersect = bvh->intersect(randomR, &randomIsect);
      if (didIntersect == true) {
        Vector3D L = randomIsect.bsdf->get_emission();
        Vector3D eval = L * Fr / pdf;
        L_out += eval;
      } 
    }
    L_out /= num_samples;
    return L_out;
  }

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  /*
  • Iterate thru all `lights`
  • `num samples = 1` if `is_delta_light()`
  • Take `num samples` # of samples per light (sampleL)
  • Intersection test (check if blocked, if not–add to L_out)
  */
  int num_samples = 0;
  Vector3D wi; Vector3D res(0.0); double dist; double pdf;
  for (auto l = scene->lights.begin(); l != scene->lights.end(); l++) {
    num_samples = (*l)->is_delta_light() ? 1 : ns_area_light;
    for (int i = 0; i < num_samples; i++) {
      Vector3D L = (*l)->sample_L(hit_p, &wi, &dist, &pdf);
      /*if (wi.z < 0) {
        continue;
      } */
      Ray randomR = Ray(hit_p, wi);
      randomR.min_t = EPS_F;
      Intersection randomIsect;
      randomR.max_t = dist - EPS_F;
      bool didIntersect = bvh->intersect(randomR, &randomIsect);
      
      if (!didIntersect) {
        Vector3D Fr = isect.bsdf->f(w_out, wi);
        Vector3D eval = L * Fr / pdf;
        L_out += eval;
      }
    }
    L_out /= num_samples;
    res += L_out;
  }

  return res;

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light

  //return Vector3D(1.0);
 return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  if(direct_hemisphere_sample) {
    return estimate_direct_lighting_hemisphere(r, isect);
  } else {
    return estimate_direct_lighting_importance(r, isect);
  }
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  const double CPDF = 0.4;
  Vector3D sample_vec; double pdf;
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();
  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out = one_bounce_radiance(r, isect);
  Vector3D Fr = isect.bsdf->sample_f(w_out, &sample_vec, &pdf); // this function gives us a random vector and calculates Fr
  sample_vec = o2w * sample_vec;
  Ray randomR = Ray(hit_p, -1 * sample_vec); // wasn't working unless we use the ray constructor??
  randomR.min_t = EPS_F;
  Intersection randomIsect;
  randomR.depth = r.depth-1;
  bool shouldStop = coin_flip(CPDF);
  bool didIntersect = bvh->intersect(randomR, &randomIsect);
  if(!shouldStop && didIntersect && randomR.depth > 1) {
    L_out+= at_least_one_bounce_radiance(r, randomIsect) * Fr / pdf / CPDF;
  }
  return L_out;

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.

  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;

  L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);//(zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect));

  // TODO (Part 3): Return the direct illumination.

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
  Vector3D estimate = Vector3D(0, 0, 0);

  for (int i = 0; i < num_samples; i++) {
    Vector2D sample = gridSampler->get_sample();
    Vector2D sample_pt = origin + sample;
    Ray ray = camera->generate_ray(sample_pt.x/sampleBuffer.w, sample_pt.y/sampleBuffer.h);
    ray.depth = this->max_ray_depth;
    estimate += est_radiance_global_illumination(ray);
    // no need to divide by pdf, since we would just be dividing by (1*1)
  }
  estimate = estimate / num_samples;

  sampleBuffer.update_pixel(estimate, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;


}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
