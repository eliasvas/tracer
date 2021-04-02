#ifndef MATERIAL_H
#define MATERIAL_H
#include "tools.h"
#include "hitable.h"

internal vec3 
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
  DIELECTRIC = 3,
}MaterialType;

typedef struct LambertianMaterial
{
  vec3 albedo;
}LambertianMaterial;

internal i32 lambertian_scatter(LambertianMaterial *m, Ray r, HitRecord *rec, vec3 *attenuation, Ray *scattered)
{
  vec3 target = vec3_add(rec->p, vec3_add(rec->normal, random_in_unit_sphere())); 
  *scattered = ray_init(rec->p, vec3_sub(target, rec->p), RAY_REFLECTION);
  *attenuation = m->albedo;
  return 1;
}


typedef struct MetalMaterial
{
  vec3 albedo;
  f32 fuzz;
}MetalMaterial;

internal vec3 vec3_reflect(vec3 v, vec3 n)
{
  return vec3_sub(v, vec3_mulf(n, vec3_dot(v,n)*2.f));
}

internal i32 metal_scatter(MetalMaterial *m, Ray r, HitRecord *rec, vec3 *attenuation, Ray *scattered)
{
  vec3 reflected = vec3_reflect(vec3_normalize(r.d),rec->normal); 
  *scattered = ray_init(rec->p, vec3_add(reflected, vec3_mulf(random_in_unit_sphere(),m->fuzz)), RAY_REFLECTION);
  *attenuation = m->albedo;
  return (vec3_dot(scattered->d, rec->normal) > 0);
}
internal i32 refract(vec3 v, vec3 n, f32 ni_over_nt, vec3* refracted)
{
  vec3 uv = vec3_normalize(v);
  f32 dt = vec3_dot(uv,n);
  f32 d = 1.f - ni_over_nt * (1.f - dt*dt);
  if (d >0)
  {
    *refracted = vec3_sub(vec3_mulf(vec3_sub(uv, vec3_mulf(n,dt)), ni_over_nt), vec3_mulf(n,sqrt(d)), RAY_REFRACTION);
    return 1;
  }
  return 0;
}


typedef struct DielectricMaterial
{
  f32 ref_index;
}DielectricMaterial;

internal i32 dielectric_scatter(DielectricMaterial *m, Ray r, HitRecord *rec, vec3 *attenuation, Ray *scattered)
{
  vec3 outward_normal;
  vec3 reflected = vec3_reflect(r.d, rec->normal);
  f32 ni_over_nt;
  *attenuation = v3(1.f,1.f,1.f);
  vec3 refracted;
  if (vec3_dot(r.d,rec->normal) > 0)
  {
    outward_normal = vec3_mulf(rec->normal, -1.f);
    ni_over_nt = m->ref_index;
  }
  else
  {
    outward_normal = rec->normal;
    ni_over_nt = 1.f / m->ref_index;
  }

  if (refract(r.d, outward_normal, ni_over_nt, &refracted))
  {
    *scattered = ray_init(rec->p, refracted, RAY_REFRACTION);
  }
  else 
  {
      *scattered = ray_init(rec->p, reflected, RAY_REFLECTION);
      return 0;
  }
  return 1;
}
typedef struct Material
{
  union
  {
    LambertianMaterial lm;
    MetalMaterial mm;
    DielectricMaterial dm;
  };
  MaterialType type;
}Material;

internal i32 scatter(Material *m,Ray r, HitRecord *rec,vec3 *attenuation,Ray* scattered)
{
  if (m->type == LAMBERTIAN)
    return lambertian_scatter(&m->lm,r, rec, attenuation, scattered);
  else if (m->type == METAL)
    return metal_scatter(&m->mm,r, rec, attenuation, scattered);
  else if (m->type == DIELECTRIC)
    return dielectric_scatter(&m->dm,r, rec, attenuation, scattered);
  return 0; 
    
}
#endif
