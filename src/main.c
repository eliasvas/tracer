#include "tools.h"
#include "ray.h"
#include "hitable.h"
#include "camera.h"
#include "material.h"
/*TODO
    -[]replace qsort!
    -[]Handedness changes for triangles?
    -[]Make ppm output with just one fwrite statement (-maximum cross platformness)
*/

internal u32 window_width = 800;
internal u32 window_height= 400;
internal u32 SAMPLES_PER_PIXEL = 4;
internal BVHNode root;
internal vec3 color(Ray r, HitableList *world, i32 depth)
{
  HitRecord rec;

  //if (hitable_list_hit(world, r, 0.01f, 30000.f, &rec))
  if (bvh_hit(root, r, 0.01f, 30000.f, &rec))
  {
    Ray scattered;
    vec3 attenuation;
    if (depth <10 && scatter(rec.m,r, &rec, &attenuation, &scattered))
      return vec3_mul(attenuation, color(scattered, world, depth+1));
    else
      return v3(0,0,0);
  }
  vec3 unit_direction = vec3_normalize(r.d);
  f32 t = 0.5f * (unit_direction.y + 1.f);
  return vec3_add(vec3_mulf(v3(1.f,1.f,1.f),1.f - t), vec3_mulf(v3(0.5, 0.4f,0.7f), t));

}
internal i32 
main(void)
{
  vec3 *framebuffer = ALLOC(sizeof(vec3) * window_width * window_height);
  seed_random((int)framebuffer); //framebufferframebuffer  is a random number
  //Camera cam = camera_lookat(v3(5,5,1),v3(0,0,-1), v3(0,1,0), 45.f, window_width / (f32)window_height, 0,1);
  Camera cam = camera_lookat(v3(0,2,0),v3(0,0,-1), v3(0,1,0), 90.f, window_width / (f32)window_height, 0,1);
  Hitable *list[8];
  list[0] = ALLOC(sizeof(Hitable));
  list[0]->type = TRIANGLE;
  list[0]->t = (Triangle){v3(0,0,-3),v3(1,0,-3),v3(0,1,-2.5)};
  list[0]->m = ALLOC(sizeof(Material));
  list[0]->m->type = LAMBERTIAN;
  list[0]->m->lm = (LambertianMaterial){v3(0.5f,0.2f,0.2f)};
  list[1] = ALLOC(sizeof(Hitable));
  list[1]->type = SPHERE;
  list[1]->s = (Sphere){v3(0,-100.5f,-1.f), 100.f};
  list[1]->m = ALLOC(sizeof(Material));
  list[1]->m->type = LAMBERTIAN;
  list[1]->m->lm = (LambertianMaterial){v3(0.3f,0.5f,0.6f)};
  list[2] = ALLOC(sizeof(Hitable));
  list[2]->type = SPHERE;
  list[2]->s = (Sphere){v3(1.f,0,-1), 0.5f};
  list[2]->m = ALLOC(sizeof(Material));
  list[2]->m->type = METAL;
  list[2]->m->mm = (MetalMaterial){v3(0.9f,0.8f,0.8f), 0.1f};
  list[3] = ALLOC(sizeof(Hitable));
  list[3]->type = SPHERE;
  list[3]->s = (Sphere){v3(-1.f,0,-1), 0.5f};
  list[3]->m = ALLOC(sizeof(Material));
  list[3]->m->type = DIELECTRIC;
  list[3]->m->dm = (DielectricMaterial){1.5f};
  list[4] = ALLOC(sizeof(Hitable));
  list[4]->type = SPHERE;
  list[4]->s = (Sphere){v3(0,0,-1), 0.5f};
  list[4]->m = ALLOC(sizeof(Material));
  list[4]->m->type = LAMBERTIAN;
  list[4]->m->lm = (LambertianMaterial){v3(0.9f,0.3f,0.3f)};


  HitableList hlist;
  hlist.list = list;
  hlist.list_size = 5;

  root = construct_bvh_tree(hlist.list, hlist.list_size, 0, 1);
  printf("tracing rays\n");
  for (i32 i = 0; i < window_width * window_height; ++i)
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
      //find the interpolated color
      col = vec3_divf(col, SAMPLES_PER_PIXEL);
      //apply gamma correction
      col = v3(sqrt(col.x),sqrt(col.y),sqrt(col.z));
      //if (col.x < 0.1f || col.y < 0.1f || col.z < 0.1f)printf("color[%i, %i] = (%f %f %f)\n", x,y, col.x, col.y,col.z);
      framebuffer[i] = v3(col.x,col.y,col.z);
  }
  printf("total intersections performed: %i\n", total_intersections);
  printf("writing to disk\n");
  ppm_save_pixels(window_width, window_height, (f32*)framebuffer);
  printf("finished\n");
}
