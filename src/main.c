#include "tools.h"
#include "ray.h"

static u32 window_width = 800;
static u32 window_height= 400;

static f32 
hit_sphere(vec3 center, f32 radius, Ray r)
{
  vec3 oc = vec3_sub(r.o, center);
  f32 a = vec3_dot(r.d, r.d);
  f32 b = 2.f * vec3_dot(oc, r.d);
  f32 c = vec3_dot(oc,oc) - radius * radius;

  f32 d = b*b - 4*a*c;
  if (d < 0)
    return -1.f;
  else 
    return (-b - sqrt(d)) / (2.f * a);
}

static vec3
color(Ray r)
{
  f32 t = hit_sphere(v3(0,0,-1), 0.5, r);
  if (t > 0.f)
  {
    vec3 N = vec3_normalize(vec3_sub(ray_point_at(r, t), v3(0,0,-1)));
    return vec3_mulf(v3(N.x + 1, N.y + 1, N.z + 1), 0.5);
  }
  vec3 unit_direction = vec3_normalize(r.d);
  t = 0.5f * (unit_direction.y + 1.f);
  return vec3_add(vec3_mulf(v3(1.f,1.f,1.f),1.f - t), vec3_mulf(v3(0.5, 0.4f,0.7f), t));
}

static int 
main(void)
{
  f32 *pixels = ALLOC(sizeof(f32) * window_width * window_height * 3);
 
  vec3 lower_left_corner = v3(-2.f, -1.f, -1.f);
  vec3 horizontal = v3(4.f,0.f,0.f);
  vec3 vertical = v3(0.f,2.f,0.f);
  vec3 origin = v3(0.f,0.f,0.f);
  for (int i = 0; i < window_width * window_height; ++i)
  {
      u32 x = i % window_width;
      u32 y = i / window_width;
      f32 u = x / (f32)window_width;
      f32 v = y / (f32)window_height;
      Ray r = ray_init(origin, vec3_add(lower_left_corner, vec3_add(vec3_mulf(horizontal, u), vec3_mulf(vertical, v))));
      vec3 col = color(r);
      pixels[3*i] = col.x;
      pixels[3*i+1] = col.y;
      pixels[3*i+2] = col.z;
  }
  ppm_save_pixels(window_width, window_height, pixels);
}
