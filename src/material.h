#ifndef MATERIAL_H
#define MATERIAL_H
#include "tools.h"
#include "hitable.h"

static vec3 
random_in_unit_sphere(void)
{
  vec3 p;
  do
  {
    p = vec3_sub(vec3_mulf(v3(random01(),random01(),random01()), 2.f), v3(1,1,1)); 
  }while (vec3_length(p) * vec3_length(p) >= 1.f);
  return p;
}
typedef enum MaterialType
{
  LAMBERTIAN = 1,
  METAL = 2,
}MaterialType;

typedef struct LambertianMaterial
{
  vec3 albedo;
}LambertianMaterial;

static i32 lambertian_scatter(LambertianMaterial *m, Ray r, HitRecord *rec, vec3 *attenuation, Ray *scattered)
{
  vec3 target = vec3_add(rec->p, vec3_add(rec->normal, random_in_unit_sphere())); 
  *scattered = ray_init(rec->p, vec3_sub(target, rec->p));
  *attenuation = m->albedo;
  return 1;
}


typedef struct MetalMaterial
{
  vec3 albedo;
  f32 fuzz;
}MetalMaterial;

static vec3 vec3_reflect(vec3 v, vec3 n)
{
  return vec3_sub(v, vec3_mulf(n, vec3_dot(v,n)*2.f));
}

static i32 metal_scatter(MetalMaterial *m, Ray r, HitRecord *rec, vec3 *attenuation, Ray *scattered)
{
  vec3 reflected = vec3_reflect(vec3_normalize(r.d),rec->normal); 
  *scattered = ray_init(rec->p, vec3_add(reflected, vec3_mulf(random_in_unit_sphere(),m->fuzz)));
  *attenuation = m->albedo;
  return (vec3_dot(scattered->d, rec->normal) > 0);
}

typedef struct Material
{
  union
  {
    LambertianMaterial lm;
    MetalMaterial mm;
  };
  MaterialType type;
}Material;

static i32 scatter(Material *m,Ray r, HitRecord *rec,vec3 *attenuation,Ray* scattered)
{
  if (m->type == LAMBERTIAN)
    return lambertian_scatter(&m->lm,r, rec, attenuation, scattered);
  else if (m->type == METAL)
    return metal_scatter(&m->mm,r, rec, attenuation, scattered);
  return 0; 
    
}
#endif
