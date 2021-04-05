#ifndef TEXTURE_H
#define TEXTURE_H
#include "tools.h"

typedef struct Perlin
{
    f32 *ranfloat;
    i32 *perm_x;
    i32 *perm_y;
    i32 *perm_z;
}Perlin;

typedef enum TextureType
{
    CONSTANT_TEXTURE = 1,
    IMAGE_TEXTURE = 2,
    NOISE_TEXTURE = 3,
    CHECKER_TEXTURE = 4,
}TextureType;

typedef struct ConstantTexture
{
    vec3 color;
}ConstantTexture;
typedef struct NoiseTexture
{
    Perlin noise;
}NoiseTexture;

typedef struct Texture Texture;
typedef struct CheckerTexture
{
    Texture *odd;
    Texture *even;
}CheckerTexture;

typedef struct ImageTexture
{
    u8 *pixels;
    i32 width;
    i32 height;
}ImageTexture;
/*
typedef struct NoiseTexture
{
    Perlin noise;
}NoiseTexture;
*/
typedef struct Texture
{
    union 
    {
        ConstantTexture ct;
        CheckerTexture checker_texture;
        NoiseTexture noise_texture;
        ImageTexture image_texture;
    };
    TextureType type;
}Texture;

internal Texture image_texture_init(char *filename)
{
    Texture tex;
    TGAInfo *image = tga_load(filename);
    tex.image_texture.width = (i32)image->width;
    tex.image_texture.height = (i32)image->height;
    tex.image_texture.pixels = image->image_data;
    tex.type = IMAGE_TEXTURE;
    return tex;
}
internal Texture constant_texture_init(vec3 color)
{
    Texture tex;
    tex.ct.color = color;
    tex.type = CONSTANT_TEXTURE;
    return tex;
}
internal Texture noise_texture_init(Perlin p)
{
    Texture tex;
    tex.type = NOISE_TEXTURE;
    tex.noise_texture.noise = p;
}
internal Texture checker_texture_init(vec3 odd_color, vec3 even_color)
{
    Texture tex;
    tex.checker_texture.odd = ALLOC(sizeof(Texture));
    tex.checker_texture.even = ALLOC(sizeof(Texture));
    *(tex.checker_texture.odd) = constant_texture_init(odd_color);
    *(tex.checker_texture.even) = constant_texture_init(even_color);
    tex.type = CHECKER_TEXTURE;
    return tex;
}

internal vec3 texture_value(Texture *tex, f32 u, f32 v, vec3 p);

internal vec3 checker_color(CheckerTexture tex, f32 u, f32 v, vec3 p)
{
    f32 sines = sin(10.f * p.x) * sin(10.f * p.y) * sin(10.f * p.z);
    if (sines < 0.f)
        return texture_value(tex.odd, u, v, p);
    else
        return texture_value(tex.even, u, v, p);
}

internal vec3 color_texture(ImageTexture tex,f32 u,f32 v,vec3 p)
{
    i32 i = u * tex.width;
    i32 j =  v * tex.height - 0.001;
    if (j < 0)j = 0;
    if (i < 0)i = 0;
    if (i > tex.width - 1)i =tex.width-1;
    if (j > tex.height - 1)j = tex.height-1;
    f32 r = ((i32)(tex.pixels[4 * i + 4 * tex.width * j]))/ 255.f;
    f32 g = ((i32)(tex.pixels[4 * i + 4 * tex.width * j + 1]))/ 255.f;
    f32 b = ((i32)(tex.pixels[4 * i + 4 * tex.width * j + 2]))/ 255.f;
    return v3(r,g,b);
}
//I have added uvp arguments because they are needed for sampling in most materials they arent even used..
internal vec3 texture_value(Texture *tex, f32 u, f32 v, vec3 p)
{
    if (tex->type == CONSTANT_TEXTURE)
        return tex->ct.color;
    if (tex->type == CHECKER_TEXTURE)
        return checker_color(tex->checker_texture, u, v, p);
    if (tex->type == NOISE_TEXTURE)
        return vec3_mulf(v3(1,1,1), noise(tex->noise_texture.noise,p));
    if (tex->type == IMAGE_TEXTURE)
        return color_texture(tex->image_texture, u, v, p);

    return v3(0,0,0);
}


internal f32* perlin_generate(void)
{
   f32 *p = ALLOC(sizeof(f32) * 256); 
   for (u32 i = 0; i < 256; ++i)
       p[i] = random01();
   return p;
}

internal void permute(i32 *p, i32 n)
{
    //this performs random swaps?!??
    for (u32 i = n-1; i > 0; --i)
    {
        i32 target = (i32)(random01() * (i + 1));
        i32 tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
}

internal i32* perlin_generate_perm(void)
{
    i32 *p = ALLOC(sizeof(i32) * 256);
    for (u32 i = 0; i < 256; ++i)
        p[i] = i;
    permute(p, 256);
    return p;
}

internal Perlin perlin;
internal f32 noise(Perlin per, vec3 p)
{
    f32 u = p.x - floor(p.x);
    f32 v = p.y - floor(p.x);
    f32 w = p.z - floor(p.x);
    i32 i = ((int)(4 * p.x)) & 255;
    i32 j = ((int)(4 * p.y)) & 255;
    i32 k = ((int)(4 * p.z)) & 255;
    printf("%i\n", per.perm_x[i] ^ per.perm_y[j] ^ per.perm_z[k]);
    return per.ranfloat[per.perm_x[i] ^ per.perm_y[j] ^ per.perm_z[k]];
}


#endif
