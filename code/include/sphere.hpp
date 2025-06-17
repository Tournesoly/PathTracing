#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D {
public:
    Sphere() {
        // unit ball at the center
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material) {
        // 
        this->center = center;
        this->radius = radius;
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f L = center - r.getOrigin();
        Vector3f dir = r.getDirection().normalized();
        float tp = Vector3f::dot(L, dir); // 投影
        float d2 = L.squaredLength() - tp * tp; // 垂径
        
        //半弦长
        float td = sqrt(radius * radius - d2);
        
        //交点距离 t，区分光源在球内还是球外
        float t;
        if (L.squaredLength() < radius * radius) {// 内
            t = tp + td;
        } else {// 外
            if (tp < 0 || d2 > radius * radius) {//不相交
                return false;
            }
            t = tp - td;
        }
        
        // 检查交点是否在有效范围内
        if (t < tmin) {
            return false;
        }
        if (t >= h.getT()) {
            return false;
        }
        
        //更新信息
        Vector3f intersectionPoint = r.getOrigin() + dir * t;
        Vector3f normal = (intersectionPoint - center) / radius;
        h.set(t, material, normal);
        return true;
    }

protected:
    Vector3f center;
    float radius;
};


#endif
