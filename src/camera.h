#ifndef CAMERA_H
#define CAMERA_H
typedef struct Camera
{
  vec3 origin;
  vec3 lower_left_corner;
  vec3 horizontal;
  vec3 vertical;
  vec3 u;
  vec3 v;
  vec3 w;
  f32 lens_radius;
}Camera;



//this function is sad :(( (fix me)
internal vec3 random_in_unit_disk(void)
{
  vec3 p;
  do
  {
    p = vec3_sub(vec3_mulf(v3(random01(),random01(),0), 2.f), v3(1,1,0)); 
  }while (vec3_length(p) * vec3_length(p) >= 1.f);
  return p;

}
internal Ray 
get_ray(Camera *c, f32 u, f32 v)
{
  vec3 rd = vec3_mulf(random_in_unit_disk(), c->lens_radius);
  vec3 offset = vec3_add(vec3_mulf(c->u,rd.x), vec3_mulf(c->v,rd.y));
  return ray_init(vec3_add(c->origin, offset), vec3_sub(vec3_sub(vec3_add(c->lower_left_corner, vec3_add(vec3_mulf(c->horizontal, u), vec3_mulf(c->vertical, v))), c->origin),offset), RAY_PRIMARY); 
}
internal Camera 
camera_lookat(vec3 lookfrom, vec3 lookat, vec3 vup, f32 vfov, f32 aspect, f32 aperture, f32 focus_dist)
{
  Camera cam;
  vec3 u,v,w;

  cam.lens_radius = aperture/2.f;
  f32 theta = vfov * PI/180.f;
  //this is so because tan = opposite/adjacent  => tan = half_height / forward(1) = half_height (take y-z view of camera and compute)
  f32 half_height = tan(theta / 2.f);
  //because aspect = ww/wh
  f32 half_width = aspect * half_height;
  cam.origin = lookfrom;
  cam.w = vec3_normalize(vec3_sub(lookfrom, lookat));
  cam.u = vec3_normalize(vec3_cross(vup, cam.w));
  cam.v = vec3_cross(cam.w,cam.u);
  //llc = origin - u * ww/2 + v * wh/2 + forward
  cam.lower_left_corner = vec3_add(vec3_add(cam.origin, vec3_mulf(cam.u, -half_width * focus_dist)),vec3_add(vec3_mulf(cam.v, -half_height * focus_dist), vec3_mulf(cam.w, -1.f * focus_dist)));
  cam.horizontal = vec3_mulf(cam.u, 2* half_width * focus_dist);
  cam.vertical = vec3_mulf(cam.v, 2* half_height * focus_dist);

  return cam;
}

#endif
