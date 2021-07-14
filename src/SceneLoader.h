#pragma once

#include <glm/glm.hpp>
#include <embree3/rtcore.h>

#include "Constants.h"
#include "Scene.h"
#include "Integrator.h"

using namespace glm;
using namespace std;

class SceneLoader {

private:

    RTCDevice _embreeDevice;

    // image
    uvec2 _imageSize = uvec2(1280, 720);
    string _outputFileName = "out.png";
    
    // camera settings
    vec3 _cameraOrigin = vec3(-1.0f, 0.0f, 0.0f);
    vec3 _cameraLookAt = vec3(0.0f, 0.0f, 0.0f);
    vec3 _cameraUp = vec3(0.0f, 1.0f, 0.0f);
    float _cameraFieldOfView = 45.0f;
    
    // geometries and materials
    vector<mat4> _sphereTransforms;
    vector<material_t> _sphereMaterials;
    vector<vec3> _rawVertices;
    vector<uvec3> _indices;
    vector<vec3> _vertices;
    vector<material_t> _triMaterials;
    material_t _curMaterial = { // default material
        vec3(0.0f),             // diffuse
        vec3(0.0f),             // specular
        1.0f,                   // shininess
        vec3(0.0f),             // emission
        vec3(0.2f, 0.2f, 0.2f), // ambient
        1.0f,                   // roughness
        "phong"                 // brdf
    };
    
    // transform
    mat4 curTransform = mat4(1.0f);
    vector<mat4> _transformStack;
    
    // lights
    vector<directionalLight_t> _directionalLights;
    vector<pointLight_t> _pointLights;
    vector<quadLight_t> _quadLights;
    vec3 _curAttenuation = vec3(1.0f, 0.0f, 0.0f);
    unordered_map<int, int> _quadLightsTable;
    
    // integrator parameters
    string _integrator;
    string _imptSampling = "hemisphere";
    string _NEE = "off";
    bool _lightstratify = false;
    bool _RR = false;
    int _maxDepth = 5;
    int _lightsamples = 1;
    int _spp = 1;
    float _gamma = 1.0f;

public:

    SceneLoader(RTCDevice embreeDevice);
    
    vec3 loadVec3(const vector<string>& arguments, size_t startIndex = 0);
    
    uvec3 loadUVec3(const vector<string>& arguments, size_t startIndex = 0);

    void loadTriangles(vec3& v1, vec3& v2, vec3& v3, material_t& material);
    
    void executeCommand(const string& command, const vector<string>& arguments);
    
    void loadSceneData(const string& filePath);
    
    Integrator* createIntegrator();
    
    void loadEmbreeTriangles(RTCScene embreeScene);
    
    void loadEmbreeSpheres(RTCScene embreeScene);
    
    RTCScene createEmbreeScene();
    
    Scene* commitSceneData();

};

void loadScene(const string& filePath, RTCDevice device, Scene** scene);
