#ifndef LIGHT_H
#define LIGHT_H

#include <Vector3f.h>
#include "object3d.hpp"

class Light {
public:
    Light() = default;

    virtual ~Light() = default;

    virtual void getIllumination(const Vector3f &p, Vector3f &dir, Vector3f &col, float &lightDist) const = 0;
    
    Vector3f position;
    Vector3f normal; // 用于nee 面光源的normal
    float area; //同 用于 nee

};


class DirectionalLight : public Light {
public:
    DirectionalLight() = delete;

    DirectionalLight(const Vector3f &d, const Vector3f &c) {
        direction = d.normalized();
        color = c;
    }

    ~DirectionalLight() override = default;

    ///@param p unsed in this function
    ///@param distanceToLight not well defined because it's not a point light
    void getIllumination(const Vector3f &p, Vector3f &dir, Vector3f &col, float &lightDist) const override {
        // the direction to the light is the opposite of the
        // direction of the directional light source
        dir = -direction;
        col = color;
        lightDist = std::numeric_limits<float>::infinity(); // 方向光距离设为无穷大
    }

private:

    Vector3f direction;
    Vector3f color;
    Vector3f normal; // 用于nee 面光源的normal
    float area; //同 用于 nee

};

class PointLight : public Light {
public:
    PointLight() = delete;

    PointLight(const Vector3f &p, const Vector3f &c) {
        position = p;
        color = c;
    }

    PointLight(const Vector3f &p, const Vector3f &c, const Vector3f &n, float &a):position(p), color(c), normal(n), area(a) {

    }

    ~PointLight() override = default;

    void getIllumination(const Vector3f &p, Vector3f &dir, Vector3f &col, float &lightDist) const override {
        // the direction to the light is the opposite of the
        // direction of the directional light source
        dir = (position - p);
        lightDist = dir.length();
        //printf("lightDist: %f\n", lightDist);
        dir = dir / dir.length();
        col = color;
    }



    Vector3f position;
    Vector3f color;
    Vector3f normal; // 用于nee 面光源的normal
    float area; //同 用于 nee

};

#endif // LIGHT_H
