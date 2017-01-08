#include "r3d.hpp"

#include <cassert>
#include <iostream>

namespace r3d {

constexpr Real EPSILON = 1e-10;

template <typename T> R3D_INLINE T max2(T a, T b) {
  return (b > a) ? (b) : (a);
}

R3D_INLINE Real rel_diff_with_floor(Real a, Real b, Real floor = EPSILON) {
  Real am = fabs(a);
  Real bm = fabs(b);
  if (am <= floor && bm <= floor)
    return 0.0;
  return fabs(b - a) / max2(am, bm);
}

R3D_INLINE bool are_close(Real a, Real b, Real tol = EPSILON,
                          Real floor = EPSILON) {
  return rel_diff_with_floor(a, b, floor) <= tol;
}

template <Int n>
R3D_INLINE bool are_close(Vector<n> a, Vector<n> b, Real tol = EPSILON,
                          Real floor = EPSILON) {
  for (Int i = 0; i < n; ++i)
    if (!are_close(a[i], b[i], tol, floor))
      return false;
  return true;
}

} // end namespace r3d

static void test_3d() {
  r3d::Few<r3d::Vector<3>, 4> verts = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto faces = r3d::faces_from_verts(verts);
  assert(r3d::are_close(faces[0].n, -r3d::normalize(r3d::vector_3(1, 1, 1))));
  assert(r3d::are_close(faces[0].d, -faces[0].n[0]));
  for (r3d::Int i = 0; i < 3; ++i) {
    auto v = r3d::vector_3(0, 0, 0);
    v[i] = 1;
    assert(r3d::are_close(faces[i + 1].n, v));
    assert(r3d::are_close(faces[i + 1].d, 0));
  }
  r3d::Polytope<3> a;
  r3d::init(a, verts);
  assert(a.nverts == 4);
  auto b = a;
  r3d::clip(b, r3d::Few<r3d::Plane<3>, 1>({{{1, 0, 0}, -0.5}}));
  assert(b.nverts == 4);
  auto volume = r3d::measure(b);
  assert(r3d::are_close(volume, r3d::cube(0.5) / 6.0));
  auto c = a;
  r3d::clip(c, r3d::Few<r3d::Plane<3>, 1>({{{-1, 0, 0}, 0.5}}));
  assert(c.nverts == 6);
  volume = r3d::measure(c);
  assert(r3d::are_close(volume, (1. / 6.) - (r3d::cube(0.5) / 6.)));
  r3d::Few<r3d::Vector<3>, 4> verts2 = {
      {1, 0, 0}, {1, 1, 0}, {0, 0, 0}, {1, 0, 1}};
  r3d::Polytope<3> intersection;
  r3d::intersect_simplices(intersection, verts, verts2);
  volume = r3d::measure(intersection);
  assert(r3d::are_close(volume, (1. / 3.) * (1. / 4.) * (1. / 2.)));
  r3d::Few<r3d::Vector<3>, 4> far_verts = {
      {0, 0, 4}, {1, 0, 4}, {0, 1, 4}, {0, 0, 5}};
  r3d::Polytope<3> null_intersection;
  r3d::intersect_simplices(null_intersection, verts, far_verts);
  assert(null_intersection.nverts == 0);
}

static void test_2d() {
  r3d::Few<r3d::Vector<2>, 3> verts = {{0, 0}, {1, 0}, {0, 1}};
  auto faces = r3d::faces_from_verts(verts);
  assert(r3d::are_close(faces[0].n, r3d::vector_2(0, 1)));
  assert(r3d::are_close(faces[0].d, 0));
  assert(r3d::are_close(faces[1].n, -r3d::normalize(r3d::vector_2(1, 1))));
  assert(r3d::are_close(faces[1].d, -faces[1].n[0]));
  assert(r3d::are_close(faces[2].n, r3d::vector_2(1, 0)));
  assert(r3d::are_close(faces[2].d, 0));
  r3d::Polytope<2> a;
  r3d::init(a, verts);
  assert(a.nverts == 3);
  auto b = a;
  r3d::clip(b, r3d::Few<r3d::Plane<2>, 1>({{{1, 0}, -0.5}}));
  assert(b.nverts == 3);
  auto area = r3d::measure(b);
  assert(r3d::are_close(area, r3d::square(0.5) / 2.0));
  auto c = a;
  r3d::clip(c, r3d::Few<r3d::Plane<2>, 1>({{{-1, 0}, 0.5}}));
  assert(c.nverts == 4);
  area = r3d::measure(c);
  assert(r3d::are_close(area, (1. / 2.) - (r3d::square(0.5) / 2.0)));
  r3d::Few<r3d::Vector<2>, 3> verts2 = {{0, 0}, {1, 0}, {1, 1}};
  r3d::Polytope<2> intersection;
  r3d::intersect_simplices(intersection, verts, verts2);
  area = r3d::measure(intersection);
  assert(r3d::are_close(area, 1. / 4.));
  r3d::Few<r3d::Vector<2>, 3> far_verts = {{0, 4}, {1, 4}, {0, 5}};
  r3d::Polytope<2> null_intersection;
  r3d::intersect_simplices(null_intersection, verts, far_verts);
  assert(null_intersection.nverts == 0);
}

int main() {
  test_2d();
  test_3d();
}
