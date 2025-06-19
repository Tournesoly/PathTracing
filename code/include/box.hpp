#ifndef BOX_HPP
#define BOX_HPP

#include <ray.hpp>
#include <vecmath.h>
#include <vector>
#include <object3d.hpp>
#include <algorithm>

// 包围体积类
class BoundingVolume {
public:
    Vector3f lowerBound; // 包围盒的下限坐标
    Vector3f upperBound; // 包围盒的上限坐标

    BoundingVolume() : lowerBound(INFINITY, INFINITY, INFINITY), upperBound(-INFINITY, -INFINITY, -INFINITY) {}
    BoundingVolume(const Vector3f& lower, const Vector3f& upper) : lowerBound(lower), upperBound(upper) {}

    // 扩展包围盒以包含新点
    void expandToInclude(const Vector3f& point) {
        lowerBound = Vector3f(
            std::min(lowerBound.x(), point.x()),
            std::min(lowerBound.y(), point.y()),
            std::min(lowerBound.z(), point.z())
        );
        upperBound = Vector3f(
            std::max(upperBound.x(), point.x()),
            std::max(upperBound.y(), point.y()),
            std::max(upperBound.z(), point.z())
        );
    }

    // 扩展包围盒以包含另一个包围盒
    void expandToInclude(const BoundingVolume& other) {
        expandToInclude(other.lowerBound);
        expandToInclude(other.upperBound);
    }

    // 检测射线与包围盒的相交
    bool intersect(const Ray& ray, float& tHit) const {
        float tMin = -INFINITY;
        float tMax = INFINITY;

        for (int dim = 0; dim < 3; ++dim) {
            float t0 = (lowerBound[dim] - ray.origin[dim]) / ray.direction[dim];
            float t1 = (upperBound[dim] - ray.origin[dim]) / ray.direction[dim];

            if (ray.direction[dim] < 0) std::swap(t0, t1); // 调整近远点

            tMin = std::max(tMin, t0);
            tMax = std::min(tMax, t1);

            if (tMax <= tMin || tMax < 0) return false; // 无交点
        }

        tHit = tMin > 0 ? tMin : tMax; // 返回最近有效交点
        return true;
    }

    // 获取所有顶点
    std::vector<Vector3f> getVertices() const {
        std::vector<Vector3f> corners;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    corners.push_back(Vector3f(
                        lowerBound.x() + i * (upperBound.x() - lowerBound.x()),
                        lowerBound.y() + j * (upperBound.y() - lowerBound.y()),
                        lowerBound.z() + k * (upperBound.z() - lowerBound.z())
                    ));
                }
            }
        }
        return corners;
    }


    // 添加质心计算方法
    Vector3f centroid() const {
        return (lowerBound + upperBound) * 0.5f;
    }
};

#endif