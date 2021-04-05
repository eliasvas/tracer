#ifndef HITABLE_H
#define HITABLE_H
#include "tools.h"
#include "ray.h"

internal u32 total_intersections = 0;
typedef struct Material Material;
typedef struct HitRecord
{
  f32 t;
  vec3 p;
  vec3 normal;
  Material *m;
}HitRecord;

typedef enum HitableType
{
  SPHERE = 1,
  TRIANGLE = 2,
  BOX = 3,
  BVH_NODE = 4,
}HitableType;



typedef struct AABB
{
    vec3 min;
    vec3 max;
}AABB;

INLINE f32 ffmin(f32 a, f32 b) { return a < b ? a : b; }
INLINE f32 ffmax(f32 a, f32 b) { return a > b ? a : b; }
internal AABB aabb_init(vec3 max, vec3 min)
{
    return (AABB){max, min};
}

internal void swap_floats(f32 *t0, f32 *t1)
{
    f32 temp = *t0;
    *t0 = *t1;
    *t1 = temp;
} 
internal i32 aabb_hit(AABB aabb, Ray r, f32 t_min, f32 t_max)
{
    total_intersections++;
    for (i32 a = 0; a < 3; ++a)
    {
        f32 invD = 1.f / r.d.elements[a];
        f32 t0 = (aabb.min.elements[a] - r.o.elements[a]) * invD;
        f32 t1 = (aabb.max.elements[a] - r.o.elements[a]) * invD;
        if (invD < 0.f)
            swap_floats(&t0, &t1);
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if (t_max <= t_min)
            return 0;
    }
    //printf("aabb hit was succesful\n");
    return 1;
}
//used to e.g calculate the bounding box of a primitive over a time interval (where there are two boxes)
internal AABB surrounding_box(AABB box0, AABB box1)
{
    vec3 small = v3(ffmin(box0.min.x, box1.min.x), ffmin(box0.min.y, box1.min.y), ffmin(box0.min.z, box1.min.z));
    vec3 big = v3(ffmax(box0.max.x, box1.max.x), ffmax(box0.max.y, box1.max.y), ffmax(box0.max.z, box1.max.z));
    //printf("L66: bounding box generated: MIN[%f %f %f], MAX[%f %f %f]\n", small.x, small.y, small.z, big.x, big.y, big.z);
    return aabb_init(small, big);
}

typedef struct Sphere
{
  vec3 center;
  f32 radius;
}Sphere;

internal i32 sphere_bounding_box(Sphere s,f32 t0, f32 t1, AABB *box)
{
    *box = aabb_init(vec3_sub(s.center, v3(s.radius,s.radius,s.radius)), vec3_add(s.center, v3(s.radius, s.radius, s.radius)));
    //printf("SPHERE: bounding box generated: MIN[%f %f %f], MAX[%f %f %f]\n", box->min.x, box->min.y, box->min.z, box->max.x, box->max.y, box->max.z);
    return 1;//1 means that the bounding box has been made, for shapes such as infinite planes this is not possible
}

internal vec3 get_sphere_center_motion_blur(Sphere s, Ray r)
{
    return vec3_add(s.center, v3(0,r.time/2, 0));
}
internal i32 
sphere_hit_motion_blur(Sphere s, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
  vec3 oc = vec3_sub(r.o, get_sphere_center_motion_blur(s, r));
  f32 a = vec3_dot(r.d, r.d);
  f32 b = vec3_dot(oc, r.d);
  f32 c = vec3_dot(oc,oc) - s.radius * s.radius;

  f32 d = b*b - a*c;
  if (d > 0)
  {
    f32 temp = (-b - sqrt(b*b - a*c))/a;
    if (temp < t_max && temp > t_min)
    {
      rec->t = temp;
      rec->p = ray_point_at(r, rec->t);
      rec->normal = vec3_normalize(vec3_divf(vec3_sub(rec->p, get_sphere_center_motion_blur(s, r)), s.radius));
      return 1;
    }
    temp = (-b + sqrt(b*b - 4*a*c))/a;
    if (temp < t_max && temp > t_min)
    {
      rec->t = temp;
      rec->p = ray_point_at(r, rec->t);
      rec->normal = vec3_normalize(vec3_divf(vec3_sub(rec->p, get_sphere_center_motion_blur(s, r)), s.radius));
      return 1;
    }
  }
  return 0;
}

internal i32 
sphere_hit(Sphere s, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
  total_intersections++;
  vec3 oc = vec3_sub(r.o, s.center);
  f32 a = vec3_dot(r.d, r.d);
  f32 b = vec3_dot(oc, r.d);
  f32 c = vec3_dot(oc,oc) - s.radius * s.radius;

  f32 d = b*b - a*c;
  if (d > 0)
  {
    f32 temp = (-b - sqrt(b*b - a*c))/a;
    if (temp < t_max && temp > t_min)
    {
      rec->t = temp;
      rec->p = ray_point_at(r, rec->t);
      rec->normal = vec3_normalize(vec3_divf(vec3_sub(rec->p, s.center), s.radius));
      return 1;
    }
    temp = (-b + sqrt(b*b - 4*a*c))/a;
    if (temp < t_max && temp > t_min)
    {
      rec->t = temp;
      rec->p = ray_point_at(r, rec->t);
      rec->normal = vec3_normalize(vec3_divf(vec3_sub(rec->p, s.center), s.radius));
      return 1;
    }
  }
  return 0;
}

//TODO: investigate
internal void get_sphere_uv(vec3 p, f32 *u, f32 *v)
{
    f32 phi = atan2(p.z, p.x);
    f32 theta = asin(p.y);
    *u = 1 - (phi + PI) / (2 * PI);
    *v = (theta + PI/2) / PI;
}
typedef struct Triangle 
{
  vec3 v0;
  vec3 v1;
  vec3 v2;
}Triangle;
internal i32 
triangle_hit(Triangle tri, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
    total_intersections++;
	f32 t, u, v;

	vec3 v0v1 = vec3_sub(tri.v1, tri.v0);
	vec3 v0v2 = vec3_sub(tri.v2, tri.v0);
    vec3 normal = vec3_cross(v0v1, v0v2);
	
	vec3 pvec = vec3_cross(r.d, v0v2);
	
	//f32 det = vec3_dot(pvec, v0v1);
	f32 det = vec3_dot(v0v1, pvec);
	f32 kEpsilon = 0.00000001;

	// if the determinant is negative the triangle is backfacing
	// if the determinant is close to 0, the ray misses the triangle
	if (fabs(det) < kEpsilon) return 0;

	f32 invDet = 1 / det;
	
	vec3 tvec = vec3_sub(r.o,tri.v0);
	u = vec3_dot(tvec, pvec) * invDet;
	
	if (u < 0 || u > 1) return 0;
    
	vec3 qvec = vec3_cross(tvec, v0v1);
	v = vec3_dot(r.d, qvec) * invDet;
	if (v < 0 || u + v > 1) return 0;

    //if (r.type == RAY_PRIMARY)
        //printf("intersection at P: %f %f %f\n", rec->p.x, rec->p.y, rec->p.z);


	t = vec3_dot(v0v2, qvec) * invDet;

	
	if (t < 0  || t > t_max || t < t_min) return 0;

	rec->p = ray_point_at(r,t);
	rec->t = t;
	rec->normal = normal;
	//rec->normal = v3(0,0,1);


	return 1;
}

internal i32 triangle_bounding_box(Triangle tri,f32 t0, f32 t1, AABB *box)
{
    vec3 min = v3(ffmin(tri.v0.x, ffmin(tri.v1.x, tri.v2.x)), 
                  ffmin(tri.v0.y, ffmin(tri.v1.y, tri.v2.y)), 
                  ffmin(tri.v0.z, ffmin(tri.v1.z, tri.v2.z)));
vec3 max = v3(ffmax(tri.v0.x, ffmax(tri.v1.x, tri.v2.x)), 
                  ffmax(tri.v0.y, ffmax(tri.v1.y, tri.v2.y)), 
                  ffmax(tri.v0.z, ffmax(tri.v1.z, tri.v2.z)));
    *box = aabb_init(min, max);
    return 1;//1 means that the bounding box has been made, for shapes such as infinite planes this is not possible
}

internal i32 box_bounding_box(AABB aabb,f32 t0, f32 t1, AABB *box)
{
    *box = aabb;
    return 1;//1 means that the bounding box has been made, for shapes such as infinite planes this is not possible
}

typedef struct Hitable Hitable;//is this legal? -I'll make it legal!!
typedef struct BVHNode
{
    Hitable *left;
    Hitable *right;
    //if we intersect we can explore the left and right nodes until we hit the leaves
    AABB box;
}BVHNode;

internal i32 bvh_hit(BVHNode node, Ray r, f32 t_min, f32 t_max, HitRecord *rec);

typedef struct Hitable
{
  union
  {
    Sphere s;
    Triangle t;
    AABB box;
    BVHNode node;
  };
  Material *m;
  HitableType type;  
}Hitable;

internal i32 
hitable_hit(Hitable *hitable, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
    if (!hitable)return 0; // -_-
    rec->m = hitable->m;
    if (hitable->type == TRIANGLE)
        return triangle_hit(hitable->t, r, t_min, t_max, rec);
    else if (hitable->type == SPHERE)
        return sphere_hit(hitable->s, r, t_min, t_max, rec);
    else if (hitable->type == BOX)
        return aabb_hit(hitable->box, r, t_min, t_max);
    else if (hitable->type == BVH_NODE)
        return bvh_hit(hitable->node,r, t_min, t_max, rec);

    return 0;
}

internal i32 hitable_bounding_box(Hitable *hitable, f32 t0, f32 t1, AABB *box)
{
    if (hitable->type == TRIANGLE)
        return triangle_bounding_box(hitable->t,t0, t1, box);
    else if (hitable->type == SPHERE)
        return sphere_bounding_box(hitable->s, t0, t1, box);
    else if (hitable->type == BOX)
        return box_bounding_box(hitable->box, t0, t1, box);
    else if (hitable->type == BVH_NODE)
        return box_bounding_box(hitable->node.box, t0, t1, box);

}

internal i32 bvh_hit(BVHNode node, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
    //for some reason aabb hit never succeeds
    //printf("min(%f, %f, %f) max(%f, %f, %f)\n", node.box.min.x, node.box.min.y, node.box.min.z, node.box.max.x, node.box.max.y, node.box.max.z);
    if (aabb_hit(node.box, r, t_min, t_max))
    {
        //printf("we got a BVH hit in bvh_hit()!!!\n");
        HitRecord left_rec, right_rec;
        i32 hit_left = hitable_hit(node.left, r, t_min, t_max, &left_rec);
        i32 hit_right = hitable_hit(node.right, r, t_min, t_max, &right_rec);
        if (hit_left && hit_right)
        {
            if (left_rec.t < right_rec.t)
                *rec = left_rec;
            else 
                *rec = right_rec;
            return 1;
        }
        else if (hit_left)
        {
            *rec = left_rec;
            return 1;
        }
        else if (hit_right)
        {
            *rec = right_rec;
            return 1;
        }
        return 0;
    }
    return 0;
}

typedef struct HitableList
{
  Hitable **list;
  u32 list_size;
}HitableList;

internal i32 
hitable_list_hit(HitableList *hl, Ray r, f32 t_min, f32 t_max, HitRecord *rec)
{
  HitRecord temp_rec;
  i32 hit_anything = 0;
  //we need f64 to not have rounding errors
  f64 closest_so_far = t_max;
  for (u32 i = 0; i < hl->list_size;++i)
  {
    if (hitable_hit(hl->list[i],r,t_min,closest_so_far, &temp_rec))
    {
      hit_anything = 1;
      closest_so_far = temp_rec.t;
      *rec = temp_rec;
    }
  }
  return hit_anything;
}

i32 box_x_compare(void *a, void *b)
{
    AABB left_box, right_box;    
    Hitable *ah = *(Hitable **)a;
    Hitable *bh = *(Hitable **)b;
    //find the bounding box for each primitive nad sort based on x axis!
    //TODO check if bounding box was made correctly (function returns 1)
    hitable_bounding_box(ah, 0, 0, &left_box);
    hitable_bounding_box(ah, 0, 0, &right_box);
    if (left_box.min.x  -  right_box.min.x < 0)
        return -1;
    else 
        return 1;

}

i32 box_y_compare(void *a, void *b)
{
    AABB left_box, right_box;    
    Hitable *ah = *(Hitable **)a;
    Hitable *bh = *(Hitable **)b;
    //find the bounding box for each primitive nad sort based on x axis!
    //TODO check if bounding box was made correctly (function returns 1)
    hitable_bounding_box(ah, 0, 0, &left_box);
    hitable_bounding_box(ah, 0, 0, &right_box);
    if (left_box.min.y - right_box.min.y < 0)
        return -1;
    else 
        return 1;

}

i32 box_z_compare(void *a, void *b)
{
    AABB left_box, right_box;    
    Hitable *ah = *(Hitable **)a;
    Hitable *bh = *(Hitable **)b;
    //find the bounding box for each primitive nad sort based on x axis!
    //TODO check if bounding box was made correctly (function returns 1)
    hitable_bounding_box(ah, 0, 0, &left_box);
    hitable_bounding_box(ah, 0, 0, &right_box);
    if (left_box.min.z - right_box.min.z < 0)
        return -1;
    else
        return 1;

}



//void qsort(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*))
internal BVHNode construct_bvh_tree(Hitable **h, i32 size, f32 time0, f32 time1) //hitable list may just have to be a Hitable after all.....
{
    BVHNode node;
    //1) randomly choose an axis
    i32 axis = (i32)(3 * random01()); 
    //2) sort based on that axis
    if (axis == 0)
    {
        qsort(h, size, sizeof(Hitable *), box_x_compare);
    }else if (axis == 1)
    {
        qsort(h, size, sizeof(Hitable *), box_y_compare);
    }else 
    {
        qsort(h, size, sizeof(Hitable *), box_z_compare);
    }
    if (size == 1)
    {
        node.left = h[0];
        node.right = h[0];
    }
    else if (size == 2)
    {
        node.left = h[0];
        node.right = h[1];
    }
    else 
    {
        //we divide in two parts and make the Hitables, current node's children
        node.left = ALLOC(sizeof(Hitable));
        node.left->node = construct_bvh_tree(h, size/2, time0, time1);
        node.left->m = NULL;
        node.left->type = BVH_NODE;
        node.right = ALLOC(sizeof(Hitable));
        node.right->node = construct_bvh_tree(h + size/2, size - size/2, time0, time1);
        node.right->m = NULL;
        node.right->type = BVH_NODE;
    }
    //3) calc the bvh node's bounding box
    AABB left_box, right_box;
    if(!hitable_bounding_box(node.left,time0, time1, &left_box)|| !hitable_bounding_box(node.right,time0, time1, &right_box))
        printf("one bounding box was not made?!??!?\n");
    node.box = surrounding_box(left_box, right_box);

    return node;
}


#endif
