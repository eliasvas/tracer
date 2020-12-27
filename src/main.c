#include "tools.h"

static u32 window_width = 640;
static u32 window_height= 360;


static int 
main(void)
{
  printf("start\n");
  f32 *pixels = ALLOC(sizeof(f32) * window_width * window_height * 3);
 
  for (int i = 0; i < window_width * window_height; ++i)
  {
      vec3 color = v3(1.0,0.5 * ((i/window_height)/(f32)window_width),0.5);
      pixels[3*i] = color.x;
      pixels[3*i+1] = color.y;
      pixels[3*i+2] = color.z;
  }
  ppm_save_pixels(window_width, window_height, pixels);
  printf("finish");
}
