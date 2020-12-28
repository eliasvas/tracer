#ifndef CAMERA_H
#define CAMERA_H
typedef struct Camera
{
  vec3 origin;
  vec3 lower_left_corner;
  vec3 horizontal;
  vec3 vertical;
}Camera;

static Ray 
get_ray(Camera *c, f32 u, f32 v)
{
  return ray_init(c->origin, vec3_add(c->lower_left_corner, vec3_add(vec3_mulf(c->horizontal, u), vec3_mulf(c->vertical, v))));
}
#endif
