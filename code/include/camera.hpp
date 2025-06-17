#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>


class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
        const Vector3f &up, int imgW, int imgH, float angle) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->angle=angle;
    }

    Ray generateRay(const Vector2f &point) override {
        // 伸缩变量
        float fx=width/(2*tan(angle/2));
        float fy=height/(2*tan(angle/2));
        // 相机空间光向，标准化
        Vector3f drc((point.x()-width/2)/fx, (height/2 - point.y())/fy, 1);
        drc=drc.normalized();
        // 相机旋转矩阵
        Matrix3f R(horizontal, -up, direction);
        // 世界空间光向
        Vector3f drw=R*drc;
        // 世界空间光线
        Ray ray(center, drw);
        return ray;
    }
    float angle;
};

#endif //CAMERA_H
