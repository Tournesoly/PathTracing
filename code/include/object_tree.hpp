// #ifndef OB_TREE_HPP
// #define OB_TREE_HPP

// #include "object3d.hpp"
// #include "box.hpp"

// // 定向包围盒树节点
// class OBNode {
// public:
//     OBNode* left;       // 左子节点
//     OBNode* right;      // 右子节点
//     BoundingVolume box; // 定向包围盒
//     int start;          // 物体索引起始索引
//     int finish;         // 物体索引结束位置

//     OBNode() : left(nullptr), right(nullptr), start(0), finish(0) {}
//     OBNode(OBNode* left, OBNode* right, const BoundingVolume& box, int start, int finish)
//         : left(left), right(right), box(box), start(start), finish(finish) {}

//     // 检查射线与包围盒相交
//     bool intersect(const Ray& ray, float& t) const {
//         return box.intersect(ray, t);
//     }
// };

// // 定向包围盒树
// class OBTree {
// public:
//     // 使用指定范围的物体构建树
//     OBTree(std::vector<Object3D*>& objects, int start, int finish)
//         : objects(objects), start(start), finish(finish) {
//         root = buildHierarchy(start, finish);
//     }

//     // 射线相交测试
//     bool intersect(const Ray& ray, Hit& hit, float tmin) const {
//         return exploreHierarchy(ray, hit, tmin, root);
//     }

//     // 获取根节点包围盒
//     BoundingVolume getEnclosingBox() const {
//         return root ? root->box : BoundingVolume();
//     }

// private:
//     OBNode* root;                 // 根节点
//     std::vector<Object3D*>& objects; // 场景中的物体列表
//     int start;                    // 物体索引起点
//     int finish;                   // 物体索引终点

//     // 构建树形结构
//     OBNode* buildHierarchy(int start, int finish) {
//         // 计算当前节点的联合包围盒
//         BoundingVolume enclosingBox;
//         for (int i = start; i < finish; ++i) {
//             BoundingVolume objBox = objects[i]->generateBoxFromObject();
//             if (i == start) {
//                 enclosingBox = objBox;
//             } else {
//                 enclosingBox.expandToInclude(objBox);
//             }
//         }

//         // 叶节点：单个物体
//         if (finish - start <= 1) {
//             return new OBNode(nullptr, nullptr, enclosingBox, start, finish);
//         }

//         // 选择分裂轴：基于包围盒的最大维度
//         Vector3f size = enclosingBox.upperBound - enclosingBox.lowerBound;
//         int splitDim = (size.x() > size.y()) ? ((size.x() > size.z()) ? 0 : 2) : ((size.y() > size.z()) ? 1 : 2);

//         // 按质心中值分割
//         int mid = (start + finish) / 2;
//         std::nth_element(objects.begin() + start, objects.begin() + mid, objects.begin() + finish,
//             [splitDim](Object3D* a, Object3D* b) {
//                 return a->generateBoxFromObject().centroid()[splitDim] < b->generateBoxFromObject().centroid()[splitDim];
//             });

//         // 递归构建子树
//         OBNode* left = buildHierarchy(start, mid);
//         OBNode* right = buildHierarchy(mid, finish);
//         return new OBNode(left, right, enclosingBox, start, finish);
//     }

//     // 递归遍历树进行相交测试
//     bool exploreHierarchy(const Ray& ray, Hit& hit, float tmin, OBNode* node) const {
//         if (!node) return false; // 空节点

//         float t;
//         if (!node->intersect(ray, t) || t < tmin) return false; // 包围盒不相交

//         bool intersectionFound = false;
//         if (!node->left && !node->right) { // 叶节点
//             if (objects[node->start]->intersect(ray, hit, tmin)) {
//                 intersectionFound = true;
//             }
//         } else { // 非叶节点
//             if (exploreHierarchy(ray, hit, tmin, node->left)) {
//                 intersectionFound = true;
//             }
//             if (exploreHierarchy(ray, hit, tmin, node->right)) {
//                 intersectionFound = true;
//             }
//         }
//         return intersectionFound;
//     }
// };

// #endif


#ifndef TREE_H
#define TREE_H

#include <ray.hpp>
#include <vecmath.h>
#include <vector>
#include <object3d.hpp>
#include <algorithm>

// 层次结构节点
// 学习原理后参考已有代码实现

class Node{
public:
    // 左右孩子
    Node *left;
    Node *right;
    // 本节点包围盒
    BoundingVolume box;
    // 物体序号区间
    int bgn;
    int end;

    Node() : left(nullptr), right(nullptr), bgn(0), end(0) {}
    Node(Node *left, Node *right, BoundingVolume box, int bgn, int end)
        : left(left), right(right), box(box), bgn(bgn), end(end) {}

    // 节点对包围盒求交
    bool intersect(const Ray &ray) const{
        float temp;
        return box.intersect(ray, temp);
    }
};

// 层次结构树
// 学习原理后参考kd树实现

class OBTree{
public:
    // 用区间内所有物体建树
    OBTree(std::vector<Object3D*> objects, int bgn, int end) : objects(std::move(objects)), bgn(bgn), end(end) {
        root = buildTree(bgn, end);
    }

    bool intersect(const Ray &ray, Hit &hit, float tmin) const {
        return intersectTree(ray, hit, tmin, root);
    }


    BoundingVolume generateBoxFromObject() const {
        return root->box;
    }

private:
    Node *root;
    // 场景内所有物品
    std::vector<Object3D*> objects;
    // 区间起点
    int bgn;
    // 区间终点
    int end;

    Node *buildTree(int bgn, int end){
        BoundingVolume box;
        for (int i = bgn; i < end; i++){
            // 区间内包围盒逐个合并
            BoundingVolume temp = objects[i]->generateBoxFromObject();
            if (i == bgn){
                box = temp;
            }else{
                // 逐次扩大本节点包围盒
                // 使用 expandToInclude 扩展当前包围盒
                box.expandToInclude(temp);
            }
        }

        // 叶节点没有孩子
        if (end - bgn == 1){
            return new Node(nullptr, nullptr, box, bgn, end);
        }else{
            // 按照极值区间大小决定分裂维度
            int axis = box.upperBound.x() - box.lowerBound.x() > box.upperBound.y() - box.lowerBound.y() ? 0 : 1;
            axis = box.upperBound.z() - box.lowerBound.z() > (axis == 0 ? box.upperBound.x() - box.lowerBound.x() : box.upperBound.y() - box.lowerBound.y()) ? 2 : axis;
            int mid = (bgn + end) / 2;
            // 中值分割
            std::nth_element(objects.begin() + bgn, objects.begin() + mid, objects.begin() + end, [axis](Object3D *a, Object3D *b) {
                return a->generateBoxFromObject().lowerBound[axis] < b->generateBoxFromObject().lowerBound[axis];
            });
            // 递归构建左右子树
            return new Node(buildTree(bgn, mid), buildTree(mid, end), box, bgn, end);
        }
    }

    bool intersectTree(const Ray &ray, Hit &hit, float tmin, Node *node) const{
        // 空节点
        if (node == nullptr){
            return false;
        }
        // 非空但一定不相交
        float temp;
        if (!node->intersect(ray)){
            return false;
        }
        // 叶节点
        bool result = false;
        if (node->left == nullptr && node->right == nullptr){
            if (objects[node->bgn]->intersect(ray, hit, tmin)){
                result = true;
            }
        }else{
        // 递归查询左右子树
            if (intersectTree(ray, hit, tmin, node->left)){
                result = true;
            }
            if (intersectTree(ray, hit, tmin, node->right)){
                result = true;
            }
        }
        return result;
    }
};

#endif
