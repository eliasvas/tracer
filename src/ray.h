#ifndef RAY_H
#define RAY_H
typedef struct Ray
{
  vec3 o;
  vec3 d;
}Ray;

static Ray
ray_init(vec3 o, vec3 d)
{
  Ray r;
  r.o = o;
  r.d = d;
  return r;
}

static vec3 
ray_point_at(Ray r, f32 t)
{
  return vec3_add(r.o, vec3_mulf(r.d, t)); 
}

#endif
