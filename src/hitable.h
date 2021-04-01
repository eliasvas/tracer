#ifndef HITABLE_H
#define HITABLE_H
#include "tools.h"
#include "ray.h"

typedef struct Material Material;
typedef struct HitRecord
{
  f32 t;
  vec3 p;
  vec3 normal;
  Material *m;
}HitRecord;

typedef enum HitableType
{
  SPHERE = 1,
  TRIANGLE = 2,
}HitableType;
typedef struct Sphere
{
  vec3 center;
  f32 radius;
}Sphere;

internal i32 
sphere_hit(Sphere s, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
  vec3 oc = vec3_sub(r.o, s.center);
  f32 a = vec3_dot(r.d, r.d);
  f32 b = vec3_dot(oc, r.d);
  f32 c = vec3_dot(oc,oc) - s.radius * s.radius;

  f32 d = b*b - a*c;
  if (d > 0)
  {
    f32 temp = (-b - sqrt(b*b - a*c))/a;
    if (temp < t_max && temp > t_min)
    {
      rec->t = temp;
      rec->p = ray_point_at(r, rec->t);
      rec->normal = vec3_normalize(vec3_divf(vec3_sub(rec->p, s.center), s.radius));
      return 1;
    }
    temp = (-b + sqrt(b*b - 4*a*c))/a;
    if (temp < t_max && temp > t_min)
    {
      rec->t = temp;
      rec->p = ray_point_at(r, rec->t);
      rec->normal = vec3_normalize(vec3_divf(vec3_sub(rec->p, s.center), s.radius));
      return 1;
    }
  }
  return 0;
}
typedef struct Triangle 
{
  vec3 v0;
  vec3 v1;
  vec3 v2;
}Triangle;

internal i32 
triangle_hit(Triangle tri, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
  //first we compute the planes (triangles) normal vector
  vec3 v0v1 = vec3_sub(tri.v1,tri.v0);
  vec3 v0v2 = vec3_sub(tri.v2,tri.v0);
  vec3 N = vec3_cross(v0v1, v0v2);
  //N = vec3_normalize(N);
  f32 area2 = vec3_length(N);

  //(1) Lets first find the point of intersection P
  f32 n_dot_ray_dir = vec3_dot(N,r.d);
  //this is in case the triangle is paralell to ray
  if (fabs(n_dot_ray_dir) < 0.0001f)
    return 0;
  f32 d = vec3_dot(N, tri.v0);

  f32 t = (vec3_dot(N,r.o) + d) / n_dot_ray_dir;
  //we don't want things behind the ray to be rendered!
  if (t < t_min || t > t_max)return 0;

  vec3 P = vec3_add(r.o, vec3_mulf(r.d, t));
  
  //now we perform inside-outside test
  vec3 C; //vector perpendiculat to triangle plane

  vec3 edge0 = vec3_sub(tri.v1, tri.v0);
  vec3 vp0 = vec3_sub(P, tri.v0);
  C = vec3_cross(edge0,vp0);
  if (vec3_dot(N,C) < 0)return 0;


  vec3 edge1 = vec3_sub(tri.v2, tri.v1);
  vec3 vp1 = vec3_sub(P, tri.v1);
  C = vec3_cross(edge1,vp1);
  if (vec3_dot(N,C) < 0)return 0;

  vec3 edge2 = vec3_sub(tri.v0, tri.v2);
  vec3 vp2 = vec3_sub(P, tri.v2);
  C = vec3_cross(edge2,vp2);
  if (vec3_dot(N,C) < 0)return 0;

  rec->t = t;
  rec->p = P;
  rec->normal = N;

  return 1;
}



typedef struct Hitable
{
  union
  {
    Sphere s;
    Triangle t;
  };
  Material *m;
  HitableType type;  
}Hitable;

internal i32 
hitable_hit(Hitable *hitable, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
    rec->m = hitable->m;
    if (hitable->type == TRIANGLE)
        return triangle_hit(hitable->t, r, t_min, t_max, rec);
    else if (hitable->type == SPHERE)
        return sphere_hit(hitable->s, r, t_min, t_max, rec);

    return 0;
}

typedef struct HitableList
{
  Hitable **list;
  u32 list_size;
}HitableList;

internal i32 
hitable_list_hit(HitableList *hl, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
  HitRecord temp_rec;
  i32 hit_anything = 0;
  //we need f64 to not have rounding errors
  f64 closest_so_far = t_max;
  for (u32 i = 0; i < hl->list_size;++i)
  {
    if (hitable_hit(hl->list[i],r,t_min,closest_so_far, &temp_rec))
    {
      hit_anything = 1;
      closest_so_far = temp_rec.t;
      *rec = temp_rec;
    }
  }
  return hit_anything;
}

#endif
