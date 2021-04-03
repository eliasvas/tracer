#ifndef RAY_H
#define RAY_H

typedef enum RayType 
{
    RAY_PRIMARY = 1,
    RAY_SHADOW,
    RAY_REFLECTION,
    RAY_REFRACTION,
}RayType;


typedef struct Ray
{
  vec3 o;
  vec3 d;
  RayType type;
  //time the ray is cast
  f32 time;
}Ray;

internal Ray
ray_init(vec3 o, vec3 d,f32 time,  RayType type)
{
  Ray r;
  r.o = o;
  r.d = d;
  r.type = type;
  r.time = time;
  return r;
}

internal vec3 
ray_point_at(Ray r, f32 t)
{
  return vec3_add(r.o, vec3_mulf(r.d, t)); 
}

#endif
