#ifndef TEXTURE_H
#define TEXTURE_H
#include "tools.h"

typedef enum TextureType
{
    CONSTANT_TEXTURE = 1,
    IMAGE_TEXTURE = 2,
    PERLIN_TEXTURE = 3,
    CHECKER_TEXTURE = 4,
}TextureType;

typedef struct ConstantTexture
{
    vec3 color;
}ConstantTexture;

typedef struct Texture Texture;
typedef struct CheckerTexture
{
    Texture *odd;
    Texture *even;
}CheckerTexture;

typedef struct Texture
{
    union 
    {
        ConstantTexture ct;
        CheckerTexture checker_texture;
    };
    TextureType type;
}Texture;

internal Texture constant_texture_init(vec3 color)
{
    Texture tex;
    tex.ct.color = color;
    tex.type = CONSTANT_TEXTURE;
    return tex;
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

//I have added uvp arguments because they are needed for sampling in most materials they arent even used..
internal vec3 texture_value(Texture *tex, f32 u, f32 v, vec3 p)
{
    if (tex->type == CONSTANT_TEXTURE)
        return tex->ct.color;
    if (tex->type == CHECKER_TEXTURE)
        return checker_color(tex->checker_texture, u, v, p);

    return v3(0,0,0);
}



#endif
