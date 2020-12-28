#include "tools.h"
#include "ray.h"
#include "hitable.h"
#include "camera.h"
#include "material.h"

static u32 window_width = 800;
static u32 window_height= 400;
static u32 SAMPLES_PER_PIXEL = 8;

static vec3 color(Ray r, HitableList *world, i32 depth)
{
  HitRecord rec;
  if (hitable_list_hit(world, r, 0.001f, 10000.f, &rec))
  {
    Ray scattered;
    vec3 attenuation;
    if (depth <5 && scatter(rec.m,r, &rec, &attenuation, &scattered))
      return vec3_mul(attenuation, color(scattered, world, depth+1));
    else
      return v3(0,0,0);
  }
  vec3 unit_direction = vec3_normalize(r.d);
  f32 t = 0.5f * (unit_direction.y + 1.f);
  return vec3_add(vec3_mulf(v3(1.f,1.f,1.f),1.f - t), vec3_mulf(v3(0.5, 0.4f,0.7f), t));

}
static i32 
main(void)
{
  f32 *pixels = ALLOC(sizeof(f32) * window_width * window_height * 3);
  seed_random();
  vec3 lower_left_corner = v3(-2.f, -1.f, -1.f);
  vec3 horizontal = v3(4.f,0.f,0.f);
  vec3 vertical = v3(0.f,2.f,0.f);
  vec3 origin = v3(0.f,0.f,0.f);
  Camera cam = (Camera){origin, lower_left_corner, horizontal, vertical};
  Hitable *list[4];
  list[0] = ALLOC(sizeof(Hitable));
  list[0]->type = SPHERE;
  list[0]->s = (Sphere){v3(0,0,-1), 0.5f};
  list[0]->m = ALLOC(sizeof(Material));
  list[0]->m->type = LAMBERTIAN;
  list[0]->m->lm = (LambertianMaterial){v3(0.9f,0.5f,0.5f)};
  list[1] = ALLOC(sizeof(Hitable));
  list[1]->type = SPHERE;
  list[1]->s = (Sphere){v3(0,-100.5f,-1.f), 100.f};
  list[1]->m = ALLOC(sizeof(Material));
  list[1]->m->type = LAMBERTIAN;
  list[1]->m->lm = (LambertianMaterial){v3(0.5f,0.5f,0.5f)};
  list[2] = ALLOC(sizeof(Hitable));
  list[2]->type = SPHERE;
  list[2]->s = (Sphere){v3(1.f,0,-1), 0.5f};
  list[2]->m = ALLOC(sizeof(Material));
  list[2]->m->type = METAL;
  list[2]->m->mm = (MetalMaterial){v3(0.2f,0.2f,0.2f), 0.1f};
  list[3] = ALLOC(sizeof(Hitable));
  list[3]->type = SPHERE;
  list[3]->s = (Sphere){v3(-1.f,0,-1), 0.5f};
  list[3]->m = ALLOC(sizeof(Material));
  list[3]->m->type = METAL;
  list[3]->m->mm = (MetalMaterial){v3(0.9f,0.8f,0.8f),0.f};


  HitableList hlist;
  hlist.list = list;
  hlist.list_size = 4;
  printf("tracing rays\n");
  for (int i = 0; i < window_width * window_height; ++i)
  {
      u32 x = i % window_width;
      u32 y = i / window_width;
      vec3 col = v3(0,0,0);
      for (i32 sample = 0; sample < SAMPLES_PER_PIXEL; ++sample)
      {
        f32 u = (x + random01()) / (f32)window_width;
        f32 v = (y + random01()) / (f32)window_height;
        Ray r = get_ray(&cam, u, v);
        vec3 p = ray_point_at(r, 2.f);
        col = vec3_add(col, color(r,&hlist, 0));
      }
      col = vec3_divf(col, SAMPLES_PER_PIXEL);
      col = v3(sqrt(col.x),sqrt(col.y),sqrt(col.z));
      pixels[3*i] = col.x;
      pixels[3*i+1] = col.y;
      pixels[3*i+2] = col.z;
  }
  printf("writing to disk\n");
  ppm_save_pixels(window_width, window_height, pixels);
  printf("finished\n");
}
