#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>


// TODO: Implement Group - add data structure to store a list of Object*
class Group : public Object3D {

public:
    std::vector<Object3D*> group_objects;

    Group() {

    }

    explicit Group (int num_objects) {

    }

    ~Group() override {

    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        bool intersected = false;
        for (auto obj : group_objects) {
            if (obj->intersect(r, h, tmin)) {
                intersected = true;
            }
        }
        return intersected;
    }

    void addObject(int index, Object3D *obj) {
        group_objects.push_back(obj);
    }

    int getGroupSize() {
        return group_objects.size();
    }

private:

};

#endif
	
