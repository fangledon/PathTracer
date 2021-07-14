#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

#include <glm/glm.hpp>
#include <embree3/rtcore.h>

using namespace glm;
using namespace std;

struct camera_t {
    vec3 origin;
    vec3 imagePlaneTopLeft;
    vec3 pixelRight;
    vec3 pixelDown;
};

struct material_t {
    vec3 diffuse;
    vec3 specular;
    float shininess;
    vec3 emission;
    vec3 ambient;
    float roughness;
    string brdf;
};

struct directionalLight_t {
    vec3 toLight;
    vec3 brightness;
};

struct pointLight_t {
    vec3 point;
    vec3 brightness;
    vec3 attenuation;
};

struct quadLight_t {
    vec3 a;
    vec3 ab;
    vec3 ac;
    vec3 nl;
    vec3 intensity;
};

class Scene {

public:

    uvec2 imageSize;
    int maxDepth;
    int lightsamples;
    int spp; // sample per pixel in monte carlo rendering
    bool lightstratify;
    bool RR;
    float gamma;
    string NEE;
    string imptSampling;
    string outputFileName;
    string integrator;
    camera_t camera;
    vector<mat3> sphereNormalTransforms;
    vector<material_t> sphereMaterials;
    vector<material_t> triMaterials;
    vector<directionalLight_t> directionalLights;
    vector<pointLight_t> pointLights;
    vector<quadLight_t> quadLights;
    unordered_map<int, int> quadLightsTable; // hitID -> index in quadlights
    RTCScene embreeScene;

    bool castRay(vec3 origin, vec3 direction,
                 vec3* hitPosition, vec3* hitNormal,
                 material_t* hitMaterial, int* hitID = nullptr) const;

    bool castOcclusionRay(vec3 origin, vec3 direction,
                          float maxDistance = numeric_limits<float>::infinity()) const;

};
