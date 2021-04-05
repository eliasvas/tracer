#include "tools.h"
#include "ray.h"
#include "hitable.h"
#include "camera.h"
#include "material.h"
#include "time.h"
/*TODO
    -[]replace qsort!
    -[]Handedness changes for triangles?
    -[]Triangles that intersect dissapear??
    -[]Make ppm output with just one fwrite statement (-maximum cross platformness)
    -[] U and V variables should be update at the hitrecord
*/

internal u32 window_width = 800;
internal u32 window_height= 400;
internal u32 SAMPLES_PER_PIXEL = 100;
internal BVHNode root;
internal Texture red;
internal Texture blue;
internal Texture green;
internal Texture white;
internal Texture matte;
internal Texture checker_texture;
internal Texture noise_texture;
internal Texture image_texture;
#define CORNELL_WIDTH 3
#define CORNELL_HEIGHT 3
#define CORNELL_DEPTH 3
#define LIGHTS_ON 1
internal vec3 color(Ray r, HitableList *world, i32 depth)
{
  HitRecord rec;
#if LIGHTS_ON
  //if (hitable_list_hit(world, r, 0.01f, 30000.f, &rec))
  if (bvh_hit(root, r, 0.01f, 30000.f, &rec))
  {
    Ray scattered;
    vec3 attenuation;
    vec3 emitted = emit(rec.m, rec.u, rec.v, rec.p);
    if (depth <50 && scatter(rec.m,r, &rec, &attenuation, &scattered))
      return vec3_add(emitted,vec3_mul(attenuation, color(scattered, world, depth+1)));
    else
      return emitted;//for light materials there is no bounce, so they are colored like so
  }
  return v3(0,0,0);
#else
  if (hitable_list_hit(world, r, 0.01f, 30000.f, &rec))
  //if (bvh_hit(root, r, 0.01f, 30000.f, &rec))
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
#endif
}

internal i32 
main(void)
{
  perlin.ranfloat = perlin_generate();
  perlin.perm_x = perlin_generate_perm();
  perlin.perm_y = perlin_generate_perm();
  perlin.perm_z = perlin_generate_perm();
  noise_texture = noise_texture_init(perlin);
  checker_texture = checker_texture_init(v3(0.9,0.9,0.9), v3(0.2,0.7,0.1));
  image_texture = image_texture_init("../texture.tga");
  red = constant_texture_init(v3(0.8,0.2,0.2));
  white = constant_texture_init(v3(2,2,2));
  blue = constant_texture_init(v3(0.3,0.3,0.9));
  matte = constant_texture_init(v3(0.9,0.9,0.9));
  green = constant_texture_init(v3(0.2,0.9,0.2));
  clock_t clk = clock();
  vec3 *framebuffer = ALLOC(sizeof(vec3) * window_width * window_height);
  seed_random((int)framebuffer); //framebufferframebuffer  is a random number
  //Camera cam = camera_lookat(v3(5,5,1),v3(0,0,-1), v3(0,1,0), 45.f, window_width / (f32)window_height, 0,1);
  Camera cam = camera_lookat(v3(0,-1.f,3),v3(0,-1,0), v3(0,1,0), 90.f, window_width / (f32)window_height, 0,1);
  Hitable *list[100];
  list[0] = ALLOC(sizeof(Hitable));
  list[0]->type = TRIANGLE;
  list[0]->t = (Triangle){v3(0,0,-30),v3(1,0,-30),v3(0,1,-20.5)};
  list[0]->m = ALLOC(sizeof(Material));
  list[0]->m->type = LAMBERTIAN;
  list[0]->m->lm = (LambertianMaterial){&red};
  list[1] = ALLOC(sizeof(Hitable));
  list[1]->type = SPHERE;
  list[1]->s = (Sphere){v3(0,-500.5f,-1.f), 5.f};
  list[1]->m = ALLOC(sizeof(Material));
  list[1]->m->type = LAMBERTIAN;
  list[1]->m->lm = (LambertianMaterial){&checker_texture};
  list[2] = ALLOC(sizeof(Hitable));
  list[2]->type = SPHERE;
  list[2]->s = (Sphere){v3(1.5f,-2,0), 0.4f};
  list[2]->m = ALLOC(sizeof(Material));
  list[2]->m->type = METAL;
  list[2]->m->mm = (MetalMaterial){v3(0.7,0.8,0.7)};
  list[3] = ALLOC(sizeof(Hitable));
  list[3]->type = SPHERE;
  list[3]->s = (Sphere){v3(-1.5f,-1,0), 0.7f};
  list[3]->m = ALLOC(sizeof(Material));
  list[3]->m->type = DIELECTRIC;
  list[3]->m->dm = (DielectricMaterial){2.5f};
  list[4] = ALLOC(sizeof(Hitable));
  list[4]->type = SPHERE;
  list[4]->s = (Sphere){v3(0,-3,0), 0.5f};
  list[4]->m = ALLOC(sizeof(Material));
  list[4]->m->type = LAMBERTIAN;
  list[4]->m->lm = (LambertianMaterial){&red};
  list[5] = ALLOC(sizeof(Hitable));
  list[5]->type = SPHERE;
  list[5]->s = (Sphere){v3(0,4,0), 2.f};
  list[5]->m = ALLOC(sizeof(Material));
  list[5]->m->type = DIFFUSE_LIGHT;
  list[5]->m->dl = (DiffuseLightMaterial){&white};

  list[6] = ALLOC(sizeof(Hitable));
  list[6]->type = SPHERE;
  list[6]->s = (Sphere){v3(0,-500.5f,-1.f), 497.f};
  list[6]->m = ALLOC(sizeof(Material));
  list[6]->m->type = LAMBERTIAN;
  list[6]->m->lm = (LambertianMaterial){&matte};

  list[7] = ALLOC(sizeof(Hitable));
  list[7]->type = SPHERE;
  list[7]->s = (Sphere){v3(0,500.5f,-1.f), 497.f};
  list[7]->m = ALLOC(sizeof(Material));
  list[7]->m->type = LAMBERTIAN;
  list[7]->m->lm = (LambertianMaterial){&green};

  list[9] = ALLOC(sizeof(Hitable));
  list[9]->type = SPHERE;
  list[9]->s = (Sphere){v3(-500.5f, 0,-1.f), 497.f};
  list[9]->m = ALLOC(sizeof(Material));
  list[9]->m->type = LAMBERTIAN;
  list[9]->m->lm = (LambertianMaterial){&blue};

  list[8] = ALLOC(sizeof(Hitable));
  list[8]->type = SPHERE;
  list[8]->s = (Sphere){v3(500.5f, 0,-1.f), 497.f};
  list[8]->m = ALLOC(sizeof(Material));
  list[8]->m->type = LAMBERTIAN;
  list[8]->m->lm = (LambertianMaterial){&red};

  list[10] = ALLOC(sizeof(Hitable));
  list[10]->type = SPHERE;
  list[10]->s = (Sphere){v3(0,0,-500.5f), 497.f};
  list[10]->m = ALLOC(sizeof(Material));
  list[10]->m->type = LAMBERTIAN;
  list[10]->m->lm = (LambertianMaterial){&green};

  list[11] = ALLOC(sizeof(Hitable));
  list[11]->type = SPHERE;
  list[11]->s = (Sphere){v3(0,0,500.5f), 497.f};
  list[11]->m = ALLOC(sizeof(Material));
  list[11]->m->type = LAMBERTIAN;
  list[11]->m->lm = (LambertianMaterial){&blue};



  /*
  list[6] = ALLOC(sizeof(Hitable));
  list[6]->type = TRIANGLE;
  list[6]->t = (Triangle){v3(CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH),v3(CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH)};
  list[6]->m = ALLOC(sizeof(Material));
  list[6]->m->type = LAMBERTIAN;
  list[6]->m->lm = (LambertianMaterial){&red};

  list[7] = ALLOC(sizeof(Hitable));
  list[7]->type = TRIANGLE;
  list[7]->t = (Triangle){v3(CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,-CORNELL_DEPTH)};
  list[7]->m = ALLOC(sizeof(Material));
  list[7]->m->type = LAMBERTIAN;
  list[7]->m->lm = (LambertianMaterial){&blue};


  list[8] = ALLOC(sizeof(Hitable));
  list[8]->type = TRIANGLE;
  list[8]->t = (Triangle){v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH)};
  list[8]->m = ALLOC(sizeof(Material));
  list[8]->m->type = LAMBERTIAN;
  list[8]->m->lm = (LambertianMaterial){&red};


  list[9] = ALLOC(sizeof(Hitable));
  list[9]->type = TRIANGLE;
  list[9]->t = (Triangle){v3(-CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH,CORNELL_HEIGHT,-CORNELL_DEPTH), v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH)};
  list[9]->m = ALLOC(sizeof(Material));
  list[9]->m->type = LAMBERTIAN;
  list[9]->m->lm = (LambertianMaterial){&blue};


  list[10] = ALLOC(sizeof(Hitable));
  list[10]->type = TRIANGLE;
  list[10]->t = (Triangle){v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH)};
  list[10]->m = ALLOC(sizeof(Material));
  list[10]->m->type = LAMBERTIAN;
  list[10]->m->lm = (LambertianMaterial){&blue};


  list[11] = ALLOC(sizeof(Hitable));
  list[11]->type = TRIANGLE;
  list[11]->t = (Triangle){v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH), v3(-CORNELL_WIDTH, -CORNELL_HEIGHT, -CORNELL_DEPTH)};
  list[11]->m = ALLOC(sizeof(Material));
  list[11]->m->type = LAMBERTIAN;
  list[11]->m->lm = (LambertianMaterial){&blue};


  list[12] = ALLOC(sizeof(Hitable));
  list[12]->type = TRIANGLE;
  list[12]->t = (Triangle){v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH),v3(CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,-CORNELL_DEPTH)};
  list[12]->m = ALLOC(sizeof(Material));
  list[12]->m->type = LAMBERTIAN;
  list[12]->m->lm = (LambertianMaterial){&red};


  list[13] = ALLOC(sizeof(Hitable));
  list[13]->type = TRIANGLE;
  list[13]->t = (Triangle){v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,-CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,-CORNELL_DEPTH), v3(-CORNELL_WIDTH, CORNELL_HEIGHT, -CORNELL_DEPTH)};
  list[13]->m = ALLOC(sizeof(Material));
  list[13]->m->type = LAMBERTIAN;
  list[13]->m->lm = (LambertianMaterial){&red};


  list[14] = ALLOC(sizeof(Hitable));
  list[14]->type = TRIANGLE;
  list[14]->t = (Triangle){v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,-CORNELL_DEPTH)};
  list[14]->m = ALLOC(sizeof(Material));
  list[14]->m->type = LAMBERTIAN;
  list[14]->m->lm = (LambertianMaterial){&blue};


  list[15] = ALLOC(sizeof(Hitable));
  list[15]->type = TRIANGLE;
  list[15]->t = (Triangle){v3(CORNELL_WIDTH,CORNELL_HEIGHT,-CORNELL_DEPTH), v3(-CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH, CORNELL_HEIGHT, -CORNELL_DEPTH)};
  list[15]->m = ALLOC(sizeof(Material));
  list[15]->m->type = LAMBERTIAN;
  list[15]->m->lm = (LambertianMaterial){&blue};

 
  list[16] = ALLOC(sizeof(Hitable));
  list[16]->type = TRIANGLE;
  list[16]->t = (Triangle){v3(CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH)};
  list[16]->m = ALLOC(sizeof(Material));
  list[16]->m->type = LAMBERTIAN;
  list[16]->m->lm = (LambertianMaterial){&red};


  list[17] = ALLOC(sizeof(Hitable));
  list[17]->type = TRIANGLE;
  list[17]->t = (Triangle){v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(-CORNELL_WIDTH,-CORNELL_HEIGHT,CORNELL_DEPTH), v3(-CORNELL_WIDTH, CORNELL_HEIGHT, CORNELL_DEPTH)};
  list[17]->m = ALLOC(sizeof(Material));
  list[17]->m->type = LAMBERTIAN;
  list[17]->m->lm = (LambertianMaterial){&red};
*/


  /*
  list[10] = ALLOC(sizeof(Hitable));
  list[10]->type = TRIANGLE;
  list[10]->t = (Triangle){v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH),v3(CORNELL_WIDTH,CORNELL_HEIGHT,CORNELL_DEPTH)};
  list[10]->m = ALLOC(sizeof(Material));
  list[10]->m->type = DIFFUSE_LIGHT;
  list[10]->m->dl = (DiffuseLightMaterial){&red};
  */




  HitableList hlist;
  hlist.list = list;
  hlist.list_size = 12;

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
      col.x = minimum(col.x, 1.f);
      col.y = minimum(col.y, 1.f);
      col.z = minimum(col.z, 1.f);
      //apply gamma correction
      col = v3(sqrt(col.x),sqrt(col.y),sqrt(col.z));
      //if (col.x < 0.1f || col.y < 0.1f || col.z < 0.1f)printf("color[%i, %i] = (%f %f %f)\n", x,y, col.x, col.y,col.z);
      framebuffer[i] = v3(col.x,col.y,col.z);
  }
  printf("total intersections performed: %i\n", total_intersections);
  printf("writing to disk\n");
  ppm_save_pixels(window_width, window_height, (f32*)framebuffer);
  clk = clock() - clk;
  f32 end_time = (f32)clk / CLOCKS_PER_SEC;
  printf("finished in %.2f seconds\n", end_time);
}
