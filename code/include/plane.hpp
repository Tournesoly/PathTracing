#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() {

    }

    Plane(const Vector3f &normal, float d, Material *m) : Object3D(m) {
        this->normal = normal;
        this->d = d;
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        if (Vector3f::dot(normal,r.getDirection())==0){
            return false;
        }
        // 解t
        float t = (d - Vector3f::dot(normal,r.getOrigin()))/(Vector3f::dot(normal,r.getDirection()));

        if (t<tmin || t>h.getT()){
            return false;
        }

        // 更新 hit 时确保法向量朝外
        Vector3f dnormal=normal;
        if (Vector3f::dot(normal,r.getDirection())>0){
            dnormal=-normal;
        }
        h.set(t, material, dnormal);
        return true;
    }
    float d;
    Vector3f normal;
    
protected:


};

#endif //PLANE_H
		

