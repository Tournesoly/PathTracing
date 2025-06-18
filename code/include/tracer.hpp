#ifndef TRACER_HPP
#define TRACER_HPP

#include "scene_parser.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include "image.hpp"
#include "camera.hpp"

#include "group.hpp"
#include "light.hpp"
#include "utils.hpp"
#include <cmath>

using namespace std;


Vector3f reflectDir(const Vector3f& I, const Vector3f& N) {
    // 标准化入射角度
    Vector3f newI = I.normalized();
    // 返回反射方向
    return (newI - 2 * Vector3f::dot(newI, N) * N).normalized();
}


// 计算折射方向
Vector3f refractDir_whitted(const Vector3f& I, const Vector3f& N, float ior) {
    Vector3f normalized_I = I.normalized();
    Vector3f normalized_N = N.normalized();

    // 判断是从外部进入还是从内部出射
    float cosThetaI = -Vector3f::dot(normalized_I, normalized_N);

    if (cosThetaI < 0) {
        // 光线从物体内部射出 -> 调整法线方向和折射率
        return refractDir_whitted(normalized_I, -normalized_N, 1.0f / ior);
    }

    // 折射率比值
    float eta = 1.0f / ior;

    // sin^2(theta_t) = eta^2 * (1 - cos^2(theta_i))
    float sin2ThetaT = eta * eta * (1.0f - cosThetaI * cosThetaI);

    // 全反射判断
    if (sin2ThetaT > 1.0f) {
        // 发生全反射，返回反射方向
        return reflectDir(I, N);
    }

    // 计算折射方向
    float cosThetaT = sqrt(1.0f - sin2ThetaT);
    Vector3f T = eta * normalized_I + (eta * cosThetaI - cosThetaT) * normalized_N;
    return T.normalized();
}


// 计算折射方向
Vector3f refractDir(const Vector3f& I, const Vector3f& N, double ior) {
    // 标准化入射方向和法线
    Vector3f incident = I.normalized();
    Vector3f normal = N.normalized();

    // 判断光线方向：从外部进入 (cosThetaI > 0) 或内部射出 (cosThetaI < 0)
    double cosThetaI = Vector3f::dot(incident, normal);
    Vector3f adjustedNormal = normal;
    double eta = ior; // 折射率比

    if (cosThetaI < 0) {
        // 光线从内部射出，调整法线方向和折射率
        adjustedNormal = -normal;
        eta = 1.0 / ior; // 反向传播时取倒数
        cosThetaI = -cosThetaI; // 修正 cos 值
    }

    // 计算 sin^2(theta_t) 使用 Snell 定律
    double sin2ThetaI = 1.0 - cosThetaI * cosThetaI;
    double sin2ThetaT = eta * eta * sin2ThetaI;

    // 全反射检查
    if (sin2ThetaT >= 1.0) {
        // 发生全反射，返回反射方向
        return (incident + 2.0 * cosThetaI * adjustedNormal).normalized();
    }

    // 计算 cos(theta_t)
    double cosThetaT = std::sqrt(1.0 - sin2ThetaT);

    // 计算折射方向
    Vector3f refractDir = eta * incident + (eta * cosThetaI - cosThetaT) * adjustedNormal;

    // 菲涅尔反射概率 (Schlick 近似)
    double R0 = ((1.0 - ior) * (1.0 - ior)) / ((1.0 + ior) * (1.0 + ior));
    double reflectance = R0 + (1.0 - R0) * std::pow(1.0 - cosThetaI, 5.0);

    // 随机决定反射或折射
    if (RandNum() < reflectance) {
        // 选择反射
        return (incident + 2.0 * cosThetaI * adjustedNormal).normalized();
    }

    // 选择折射
    return refractDir.normalized();
}

// 漫反射方向生成函数
Vector3f diffuseDir(const Vector3f &normal) {
    // 生成两个随机数 [0, 1)
    double r1 = 2.0 * M_PI * RandNum(); // 角度 [0, 2π)
    double r2 = RandNum();              // 随机值 [0, 1)

    // 使用平方根变换生成余弦加权分布
    double r2s = std::sqrt(r2);

    // 计算球面坐标
    double x = std::cos(r1) * r2s; // x 坐标
    double y = std::sin(r1) * r2s; // y 坐标
    double z = std::sqrt(1.0 - r2); // z 坐标，余弦加权

    // 创建临时向量
    Vector3f dir(x, y, z);

    // 构造局部坐标系
    Vector3f w = normal; // 法线作为 z 轴
    Vector3f u = Vector3f::cross(std::abs(w.x()) > 0.1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0), w).normalized(); // u 轴
    Vector3f v = Vector3f::cross(w, u); // v 轴

    // 将方向转换到世界坐标系
    Vector3f worldDir = u * dir.x() + v * dir.y() + w * dir.z();

    // 确保方向在法线半球内并归一化
    if (Vector3f::dot(worldDir,normal) < 0) {
        worldDir = -worldDir; // 如果方向在下半球，反转
    }
    return worldDir.normalized();
}


// 生成 Glossy 反射方向 (Cook-Torrance 风格)
Vector3f glossyDir(const Vector3f& In, const Vector3f& normal, float shininess) {
    // 标准化输入向量
    Vector3f incident = In.normalized();
    Vector3f N = normal.normalized();

    // 将 shininess 转换为粗糙度 (roughness)，值越小越光滑
    float roughness = 1.0f / (shininess + 1.0f); // shininess = 0 -> roughness = 1 (扩散), shininess 大 -> roughness 趋近 0 (镜面)

    // GGX 采样微表面法向量 H
    float u1 = RandNum();
    float u2 = RandNum();
    float theta = std::atan(roughness * std::sqrt(u1) / std::sqrt(1.0f - u1)); // GGX 分布角度
    float phi = 2.0 * M_PI * u2; // 方位角

    float cosTheta = std::cos(theta);
    float sinTheta = std::sin(theta);
    float x = sinTheta * std::cos(phi);
    float y = sinTheta * std::sin(phi);
    float z = cosTheta;

    // 局部半向量
    Vector3f H_local(x, y, z);

    // 构造局部坐标系
    Vector3f w = N; // 法线作为 z 轴
    Vector3f u = Vector3f::cross(std::abs(w.x()) > 0.1f ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0), w).normalized();
    Vector3f v = Vector3f::cross(w, u).normalized();

    // 将 H 转换到世界坐标系
    Vector3f H = u * H_local.x() + v * H_local.y() + w * H_local.z();

    // 计算反射方向 R = I - 2 * (I · H) * H
    double dotIH = Vector3f::dot(incident, H);
    Vector3f reflectDir = incident - 2.0 * dotIH * H;

    // 确保方向在法线半球内
    if (Vector3f::dot(reflectDir, N) < 0) {
        reflectDir = -reflectDir;
    }

    return reflectDir.normalized();
}

Vector3f glossy_BRDF(const Vector3f & normal, const Vector3f &viewDir, const Vector3f &lightDir, const Vector3f &baseReflectance, float glossiness){
    Vector3f halfVector = (viewDir + lightDir).normalized();

    // 1. Fresnel 项：Schlick 近似
    float cosThetaH = std::max(0.0f, Vector3f::dot(viewDir, halfVector));
    Vector3f fresnelTerm = baseReflectance + (Vector3f(1.0f, 1.0f, 1.0f) - baseReflectance) * std::pow(std::max(0.0f, 1.0f - cosThetaH), 5.0f);

    // 2. 几何项：Smith 模型 (Schlick-GGX)
    float roughness = 1.0f - glossiness; // 粗糙度 = 1 - 光泽度
    float k = (roughness * roughness) / 2.0f;
    float cosThetaV = std::max(0.0f, Vector3f::dot(normal, viewDir));
    float cosThetaL = std::max(0.0f, Vector3f::dot(normal, lightDir));
    float geomView = cosThetaV / (cosThetaV * (1.0f - k) + k + 1e-6f);
    float geomLight = cosThetaL / (cosThetaL * (1.0f - k) + k + 1e-6f);
    float geometryTerm = geomView * geomLight;

    // 3. 分布项：GGX 微表面分布
    float cosThetaN = std::max(0.0f, Vector3f::dot(normal, halfVector));
    float roughnessSquared = roughness * roughness;
    float cosThetaNSquared = cosThetaN * cosThetaN;
    float denominator = (cosThetaNSquared * (roughnessSquared - 1.0f) + 1.0f) * (cosThetaNSquared * (roughnessSquared - 1.0f) + 1.0f) + 1e-6f;
    float distributionTerm = roughnessSquared / (M_PI * denominator);

    // 归一化并计算最终 BRDF
    float normalization = 4.0f * cosThetaV * cosThetaL + 1e-7f;
    Vector3f specular = fresnelTerm * geometryTerm * distributionTerm / normalization;

    return specular;

}

class RayCastTracer {
public:
    RayCastTracer(const SceneParser* sceneParser, const string& outputFile) : _sceneParser(sceneParser), _outputFile(outputFile){};

    void render() {
        Camera *camera = _sceneParser->getCamera();
        Image img(camera->getWidth(), camera->getHeight());
        // 循环屏幕空间的像素
        for ( int x = 0; x < camera->getWidth() ; ++x) {
            for ( int y = 0; y < camera->getHeight() ; ++y) {
                // 计算当前像素(x,y)处相机出射光线camRay
                Ray camRay = camera->generateRay(Vector2f(x, y)) ;
                Group* baseGroup = _sceneParser->getGroup() ;
                Hit hit ;
                // 判断camRay是否和场景有交点，并返回最近交点的数据，存储在hit 中
                bool isIntersect = baseGroup->intersect (camRay, hit , 0) ;
                if ( isIntersect ) {
                    Vector3f finalColor = Vector3f ::ZERO;
                    // 找到交点之后，累加来自所有光源的光强影响
                    for ( int li = 0; li < _sceneParser->getNumLights() ; ++li ) {
                        Light* light = _sceneParser->getLight( li );
                        Vector3f L, lightColor ;
                        float fuck_dist ;
                        // 获得光照强度
                        light->getIllumination(camRay. pointAtParameter( hit .getT() ) , L, lightColor, fuck_dist) ;
                        // 计算局部光强
                        finalColor += hit.getMaterial()->Shade(camRay, hit , L, lightColor) ;
                    }
                    img.SetPixel(x, y, finalColor);
                } else {
                // 不存在交点，返回背景色
                    img.SetPixel (x, y, _sceneParser->getBackgroundColor() ) ;
                }
            }
        }
        img.SaveBMP(_outputFile.c_str());

    }

private:
    const SceneParser* _sceneParser;
    const string &_outputFile;

};

// ========================
// Whitted 风格光线追踪器
// ========================
class WhittedStyleTracer {
public:
    WhittedStyleTracer(const SceneParser* sceneParser, int maxDepth, const string& outputFile): _sceneParser(sceneParser), _maxDepth(maxDepth), _outputFile(outputFile){};

    Vector3f trace(const Ray& ray, int depth){
        if (depth <= 0) {
            return Vector3f::ZERO; // 最大递归深度限制
        }

        Hit hit;
        if (!_sceneParser->getGroup()->intersect(ray, hit, 0.001f)) {
            return _sceneParser->getBackgroundColor(); // 没有交点，返回背景色
        }

        Material* material = hit.getMaterial();
        Vector3f hitPoint = ray.pointAtParameter(hit.getT());
        Vector3f normal = hit.getNormal().normalized();
        Vector3f finalColor = Vector3f::ZERO;

        if (material->type_d_rl_rr.x()){
            //漫反射，累加来自所有光源的光强影响
            for ( int li = 0; li < _sceneParser->getNumLights() ; ++li ) {
                // 点光源
                Light* light = _sceneParser->getLight( li );
                Vector3f lightDir, lightColor ;
                float lightDist;
                // 获得光照信息
                light->getIllumination(hitPoint, lightDir, lightColor, lightDist);
                Ray shadowRay(hitPoint + normal * 0.001, lightDir);
                Hit shadowHit;
                if ( _sceneParser->getGroup()->intersect(shadowRay, shadowHit, 0.001f) && shadowHit.getT() < lightDist ){
                    printf("shadowHit and lightDist: %f %f\n", shadowHit.getT(), lightDist);
                    // 存在阴影
                }else {
                    // 计算局部光强 phong 模型
                    finalColor += hit.getMaterial()->Shade(ray, hit , lightDir, lightColor) ;
                }

            }
        } else if (material->type_d_rl_rr.y()){
            //反射
            Vector3f refl_ray = reflectDir(ray.getDirection(),normal);
            Ray reflRay(hitPoint + normal * 0.001, refl_ray);
            finalColor += trace(reflRay, depth - 1);
        } else if (material->type_d_rl_rr.z()){
            //折射
            Vector3f refr_ray = refractDir_whitted(ray.getDirection(),normal, material->refr_rate);
            Ray refrRay(hitPoint - normal * 0.001, refr_ray);
            // 有可能全反射
            if (Vector3f::dot(refr_ray, normal) > 0){
                refrRay.setOrigin(hitPoint + normal * 0.001);
            }
            finalColor += trace(refrRay, depth - 1);
        }

        return finalColor;


    }
    void render(){
        Camera* camera = _sceneParser->getCamera();
        int width = camera->getWidth();
        int height = camera->getHeight();
        Image image = Image(width, height);

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Vector2f point(i, j);
                Ray ray = camera->generateRay(point);
                Vector3f color = trace(ray, _maxDepth);
                image.SetPixel(i, j, color);
            }
        }

        image.SaveBMP(_outputFile.c_str());
        printf("WhittedTracing: Finished rendering %s\n", _outputFile.c_str());
    }

private:
    const SceneParser* _sceneParser;
    int _maxDepth;
    string _outputFile;
};

// ========================
// 简单路径追踪器（基础 Path Tracer）
// ========================
class SimplePathTracer {
public:
    SimplePathTracer(const SceneParser* sceneParser,int SPP, const string& outputFile): _sceneParser(sceneParser), _SPP(SPP), _outputFile(outputFile){};

    Vector3f trace(Ray ray) {
        Vector3f color = Vector3f::ZERO;
        Vector3f throughput = Vector3f(1, 1, 1);

        while (true) {
            Hit hit;
            if (!_sceneParser->getGroup()->intersect(ray, hit, 0.001f)) {
                color += _sceneParser->getBackgroundColor() * throughput;
                break;
            }

            Material* material = hit.getMaterial();
            Vector3f hitPoint = ray.pointAtParameter(hit.getT());
            Vector3f normal = hit.getNormal().normalized();
            Vector3f input = ray.direction.normalized();
            Vector3f emitted = material->emissionColor;
            Vector3f diffuse = material->diffuseColor;

            color += diffuse * emitted * throughput;
            throughput = throughput * diffuse;

            ray.origin = hitPoint;
            double rate = RandNum();
            Vector3f type = material->type_d_rl_rr;
            Vector3f newdir;
            if (rate <= type.x()) { // 漫反射
                ray.direction = diffuseDir(normal);
            } else if (rate <= type.x() + type.y()) { // 反射
                ray.direction = reflectDir(ray.direction, normal);
            } else { // 折射
                ray.direction = refractDir(ray.direction, normal, material->refr_rate);
            }

            if (type.x() != 0) {
                throughput *= std::abs(Vector3f::dot(ray.direction, normal));
            }

            // 俄罗斯轮盘赌
            float q = 0.1f; // 固定终止概率，可动态调整
            double rouletteRand = RandNum(); // 新生成的随机数
            if (rouletteRand < q) { // 以概率 q 终止
                break;
            }
            throughput = throughput / (1 - q); // 补偿继续概率
        }

        return color;
    }

    void render(){
        Camera* camera = _sceneParser->getCamera();
        int width = camera->getWidth();
        int height = camera->getHeight();
        Image image = Image(width, height);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Vector3f finalColor(0, 0, 0);
                for (int s = 0; s < _SPP; s++) {
                    float u = (x + RandNum());
                    float v = (y + RandNum());
                    Ray ray = camera->generateRay(Vector2f(u, v));
                    finalColor += trace(ray);
                }
                finalColor[0] /= _SPP;
                finalColor[1] /= _SPP;
                finalColor[2] /= _SPP;
                image.SetPixel(x, y, finalColor);
            }
        }

        image.SaveBMP(_outputFile.c_str());
        printf("SimplePathTracing: Finished rendering %s\n", _outputFile.c_str());
    }

private:
    const SceneParser* _sceneParser;
    int _SPP = 100;
    std::string _outputFile;

};

// ========================
// 带下一事件估计（NEE）的路径追踪器
// ========================
class NEEedPathTracer {
public:
    NEEedPathTracer(const SceneParser* sceneParser,int SPP, const string& outputFile): _sceneParser(sceneParser), _SPP(SPP), _outputFile(outputFile){};

    Vector3f trace(Ray ray) {
        Vector3f color = Vector3f::ZERO;
        Vector3f throughput = Vector3f(1, 1, 1);

        while (true) {
            Hit hit;
            if (!_sceneParser->getGroup()->intersect(ray, hit, 0.001f)) {
                color += _sceneParser->getBackgroundColor() * throughput;
                break;
            }

            Material* material = hit.getMaterial();
            Vector3f hitPoint = ray.pointAtParameter(hit.getT());
            Vector3f normal = hit.getNormal().normalized();
            Vector3f input = ray.direction.normalized();
            Vector3f emitted = material->emissionColor;
            Vector3f diffuse = material->diffuseColor;

            // nee
            if (material->type_d_rl_rr.x() != 0) {
                if (emitted.x() == 0 && emitted.y() == 0 && emitted.z() == 0) {
                    int numLights = _sceneParser->getNumEmissions();
                    if (numLights > 0) {
                        int lightIdx = std::floor(RandNum() * numLights);
                        PointLight* light = _sceneParser->getLightFromEmission(lightIdx);

                        Vector3f lightDir, lightColor;
                        float lightDist;
                        light->getIllumination(hitPoint, lightDir, lightColor, lightDist);

                        // 修正影子光线方向
                        Vector3f toLight = lightDir * lightDist;
                        Ray shadowRay(hitPoint, lightDir);
                        Hit shadowHit;
                        bool occluded = _sceneParser->getGroup()->intersect(shadowRay, shadowHit, 0.005f);

                        if (!occluded || shadowHit.getT() >= lightDist - 0.05f) {
                            float cosTheta = std::max(0.0f, Vector3f::dot(normal, lightDir));
                            if (cosTheta > 0) {
                                float lightCosTheta = std::max(0.0f, Vector3f::dot(-lightDir, light->normal));
                                if (lightCosTheta > 0) {
                                    float geoTerm = cosTheta * lightCosTheta / (lightDist * lightDist + 1e-6f);
                                    float lightPdf = (light->area > 0) ? 1.0f / (numLights * light->area + 1e-6f) : 1.0f / (numLights * 1.0f + 1e-6f);
                                    color += lightColor * diffuse * geoTerm * throughput / lightPdf / M_PI;
                                }
                            }
                        }
                        delete light;
                    }
                } else {
                    // 自发光表面贡献
                    color += diffuse * emitted * throughput / M_PI;
                    // 漫反射时乘以 1/π 归一化 BRDF
                }
            }

            double rate = RandNum();
            Vector3f type = material->type_d_rl_rr;
            Vector3f newdir;
            if (material->shininess){// 处理 glossy 材质
                // newdir = glossyDir(ray.direction, normal, material->shininess);
            } else if (rate <= type.x()) { // 漫反射
                newdir = diffuseDir(normal);
            } else if (rate <= type.x() + type.y()) { // 反射
                newdir = reflectDir(ray.direction, normal);
            } else { // 折射
                newdir = refractDir(ray.direction, normal, material->refr_rate);
            }

            ray.direction = newdir;

            if (type.x() != 0) {
                throughput *= std::abs(Vector3f::dot(ray.direction, normal));
            }

            throughput = throughput * diffuse;
            ray.origin = hitPoint;

            // 俄罗斯轮盘赌
            float q = 0.1f;
            double rouletteRand = RandNum();
            if (rouletteRand < q) {
                break;
            }
            throughput = throughput / (1 - q);
        }
        return color;
    }
    void render(){
        Camera* camera = _sceneParser->getCamera();
        int width = camera->getWidth();
        int height = camera->getHeight();
        Image image = Image(width, height);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Vector3f finalColor(0, 0, 0);
                for (int s = 0; s < _SPP; s++) {
                    float u = (x + RandNum());
                    float v = (y + RandNum());
                    Ray ray = camera->generateRay(Vector2f(u, v));
                    finalColor += trace(ray);
                }
                finalColor[0] /= _SPP;
                finalColor[1] /= _SPP;
                finalColor[2] /= _SPP;
                image.SetPixel(x, y, finalColor);
            }
        }

        image.SaveBMP(_outputFile.c_str());
        printf("SimplePathTracing: Finished rendering %s\n", _outputFile.c_str());
    }

private:
    const SceneParser* _sceneParser;
    int _SPP = 100;
    std::string _outputFile;

};

// ========================
// 实现 glossy 材质的路径追踪器
// ========================
class GlossyPathTracer {
public:
    GlossyPathTracer(const SceneParser* sceneParser, int SPP, const string& outputFile): _sceneParser(sceneParser), _SPP(SPP), _outputFile(outputFile){};

    Vector3f trace(Ray ray) {
        Vector3f color = Vector3f::ZERO;
        Vector3f throughput = Vector3f(1, 1, 1);

        while (true) {
            Hit hit;
            if (!_sceneParser->getGroup()->intersect(ray, hit, 0.001f)) {
                color += _sceneParser->getBackgroundColor() * throughput;
                break;
            }

            Material* material = hit.getMaterial();
            Vector3f hitPoint = ray.pointAtParameter(hit.getT());
            Vector3f normal = hit.getNormal().normalized();
            Vector3f input = ray.direction.normalized();
            Vector3f emitted = material->emissionColor;
            Vector3f diffuse = material->diffuseColor;

            // nee
            if (material->type_d_rl_rr.x() != 0) {
                if (emitted.x() == 0 && emitted.y() == 0 && emitted.z() == 0) {
                    int numLights = _sceneParser->getNumEmissions();
                    if (numLights > 0) {
                        int lightIdx = std::floor(RandNum() * numLights);
                        PointLight* light = _sceneParser->getLightFromEmission(lightIdx);

                        Vector3f lightDir, lightColor;
                        float lightDist;
                        light->getIllumination(hitPoint, lightDir, lightColor, lightDist);

                        // 修正影子光线方向
                        Vector3f toLight = lightDir * lightDist;
                        Ray shadowRay(hitPoint, lightDir);
                        Hit shadowHit;
                        bool occluded = _sceneParser->getGroup()->intersect(shadowRay, shadowHit, 0.005f);

                        if (!occluded || shadowHit.getT() >= lightDist - 0.05f) {
                            float cosTheta = std::max(0.0f, Vector3f::dot(normal, lightDir));
                            if (cosTheta > 0) {
                                float lightCosTheta = std::max(0.0f, Vector3f::dot(-lightDir, light->normal));
                                if (lightCosTheta > 0) {
                                    float geoTerm = cosTheta * lightCosTheta / (lightDist * lightDist + 1e-6f);
                                    float lightPdf = (light->area > 0) ? 1.0f / (numLights * light->area + 1e-6f) : 1.0f / (numLights * 1.0f + 1e-6f);
                                    color += lightColor * diffuse * geoTerm * throughput / lightPdf / M_PI;
                                }
                            }
                        }
                        delete light;
                    }
                } else {
                    // 自发光表面贡献
                    color += diffuse * emitted * throughput / M_PI;
                    // 漫反射时乘以 1/π 归一化 BRDF
                }
            }

            double rate = RandNum();
            Vector3f type = material->type_d_rl_rr;
            Vector3f newdir;
            if (material->shininess){// 处理 glossy 材质
                newdir = glossyDir(ray.direction, normal, material->shininess);
                Vector3f specular = glossy_BRDF(normal, -input, newdir, diffuse, material->shininess);
                throughput = throughput * specular; // 更新通量
            } else if (rate <= type.x()) { // 漫反射
                newdir = diffuseDir(normal);
            } else if (rate <= type.x() + type.y()) { // 反射
                newdir = reflectDir(ray.direction, normal);
            } else { // 折射
                newdir = refractDir(ray.direction, normal, material->refr_rate);
            }

            ray.direction = newdir;

            if (type.x() != 0) {
                throughput *= std::abs(Vector3f::dot(ray.direction, normal));
            }

            throughput = throughput * diffuse;
            ray.origin = hitPoint;

            // 俄罗斯轮盘赌
            float q = 0.1f; // 固定终止概率，可动态调整
            double rouletteRand = RandNum(); // 新生成的随机数
            if (rouletteRand < q) { // 以概率 q 终止
                break;
            }
            throughput = throughput / (1 - q); // 补偿继续概率
        }
        return color;

    }

    void render(){
        Camera* camera = _sceneParser->getCamera();
        int width = camera->getWidth();
        int height = camera->getHeight();
        Image image = Image(width, height);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Vector3f finalColor(0, 0, 0);
                for (int s = 0; s < _SPP; s++) {
                    float u = (x + RandNum());
                    float v = (y + RandNum());
                    Ray ray = camera->generateRay(Vector2f(u, v));
                    finalColor += trace(ray);
                }
                finalColor[0] /= _SPP;
                finalColor[1] /= _SPP;
                finalColor[2] /= _SPP;
                image.SetPixel(x, y, finalColor);
            }
        }

        image.SaveBMP(_outputFile.c_str());
        printf("SimplePathTracing: Finished rendering %s\n", _outputFile.c_str());
    }

private:
    const SceneParser* _sceneParser;
    int _SPP = 100;
    std::string _outputFile;

};
#endif // TRACER_HPP