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

static i32 
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

typedef struct Hitable
{
  union
  {
    Sphere s;
  };
  Material *m;
  HitableType type;  
}Hitable;

static i32 
hitable_hit(Hitable *hitable, Ray r, f32 t_min, f32 t_max, HitRecord * rec)
{
  rec->m = hitable->m;
  if (hitable->type == SPHERE)
    return sphere_hit(hitable->s, r, t_min, t_max, rec);
  return 0;
}

typedef struct HitableList
{
  Hitable **list;
  u32 list_size;
}HitableList;

static i32 
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
