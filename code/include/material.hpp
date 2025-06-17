#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>




// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
public:

    explicit Material(const Vector3f &d_color, const Vector3f &s_color = Vector3f::ZERO, float s = 0, const Vector3f &type_d_rl_rr = Vector3f(1,0,0), float refr_rate = 0) :
            diffuseColor(d_color), specularColor(s_color), shininess(s), type_d_rl_rr(type_d_rl_rr), refr_rate(refr_rate) {

    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;

        // diffuse
        Vector3f n=hit.getNormal();
        Vector3f l=dirToLight;

        float ln=Vector3f::dot(l,n);
        if (ln<0) {
            ln=0;
        }

        Vector3f diffuse = diffuseColor * ln;

        // specular
        Vector3f r=2*(Vector3f::dot(n,l))*n-l;
        Vector3f v=-ray.getDirection();

        float vr=Vector3f::dot(v,r);
        if (vr<0){vr=0;}

        Vector3f specular = specularColor * pow(vr, shininess);

        // 计算光照
        shaded=lightColor*(diffuse+specular);
        return shaded;
    }


    Vector3f diffuseColor; // 漫反射颜色
    Vector3f specularColor; // 镜面反射颜色
    float shininess; // 高光指数

    // whitted style 漫反射 反射 折射
    Vector3f type_d_rl_rr;

    float refr_rate;


    
    
};


#endif // MATERIAL_H
