#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <stdexcept>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "SceneLoader.h"

SceneLoader::SceneLoader(RTCDevice embreeDevice)
    : _embreeDevice(embreeDevice)
{
}

vec3 SceneLoader::loadVec3(const vector<string>& arguments, size_t startIndex)
{
    return vec3(
        stof(arguments[startIndex]),
        stof(arguments[startIndex + 1]),
        stof(arguments[startIndex + 2]));
}

uvec3 SceneLoader::loadUVec3(const vector<string>& arguments, size_t startIndex)
{
    return uvec3(
        stoi(arguments[startIndex]),
        stoi(arguments[startIndex + 1]),
        stoi(arguments[startIndex + 2]));
}

void SceneLoader::loadTriangles(vec3& v1, vec3& v2, vec3& v3, material_t& material) 
{
    _indices.push_back(uvec3(_vertices.size(), _vertices.size() + 1, _vertices.size() + 2));

    _vertices.push_back(vec3(curTransform * vec4(v1, 1.0f)));
    _vertices.push_back(vec3(curTransform * vec4(v2, 1.0f)));
    _vertices.push_back(vec3(curTransform * vec4(v3, 1.0f)));

    _triMaterials.push_back(material);
}

void SceneLoader::executeCommand(const string& command, const vector<string>& arguments)
{
    if (command == "size") {

        _imageSize = uvec2(stoi(arguments[0]), stoi(arguments[1]));

    } else if (command == "integrator") {
        
        _integrator = arguments[0];
        
    } else if (command == "spp") {
        
        _spp = stoi(arguments[0]);
        
    } else if (command == "nexteventestimation") {
        
        _NEE = arguments[0];
        
    } else if (command == "importancesampling") {
        
        _imptSampling = arguments[0];
        
    } else if (command == "russianroulette") {
        
        _RR = (arguments[0] == "on") ? true : false;
        
    } else if (command == "maxdepth") {

        _maxDepth = stoi(arguments[0]);
        if (_maxDepth == -1) {
            _maxDepth = numeric_limits<int>::max();
        }

    } else if (command == "gamma") {
        
        _gamma = stof(arguments[0]);
        
    } else if (command == "output") {

        _outputFileName = arguments[0];

    } else if (command == "camera") {

        _cameraOrigin = loadVec3(arguments, 0);
        _cameraLookAt = loadVec3(arguments, 3);
        _cameraUp = loadVec3(arguments, 6);
        _cameraFieldOfView = stof(arguments[9]);

    } else if (command == "sphere") {

        vec3 center = loadVec3(arguments, 0);
        float radius = stof(arguments[3]);

        mat4 transform = curTransform;
        transform = translate(transform, center);
        transform = scale(transform, vec3(radius));
        _sphereTransforms.push_back(transform);

        _sphereMaterials.push_back(_curMaterial);

    } else if (command == "maxverts") {

        // ignore since we are using vector

    } else if (command == "vertex") {

        _rawVertices.push_back(loadVec3(arguments));

    } else if (command == "tri") {

        uvec3 idx = loadUVec3(arguments);

        loadTriangles(_rawVertices[idx[0]], 
                      _rawVertices[idx[1]], 
                      _rawVertices[idx[2]],
                      _curMaterial);

    } else if (command == "translate") {

        vec3 translation = loadVec3(arguments);
        curTransform = translate(curTransform, translation);

    } else if (command == "rotate") {

        vec3 axis = loadVec3(arguments, 0);
        float radians = stof(arguments[3]) * PI / 180.0f;
        curTransform = rotate(curTransform, radians, axis);

    } else if (command == "scale") {

        vec3 scaleVec = loadVec3(arguments);
        curTransform = scale(curTransform, scaleVec);

    } else if (command == "pushTransform") {

        _transformStack.push_back(curTransform);

    } else if (command == "popTransform") {

        curTransform = _transformStack.back();
        _transformStack.pop_back();

    } else if (command == "directional") {

        directionalLight_t light;
        light.toLight = normalize(loadVec3(arguments, 0));
        light.brightness = loadVec3(arguments, 3);
        _directionalLights.push_back(light);

    } else if (command == "point") {

        pointLight_t light;
        light.point = loadVec3(arguments, 0);
        light.brightness = loadVec3(arguments, 3);
        light.attenuation = _curAttenuation;
        _pointLights.push_back(light);

    } else if (command == "quadLight") {
        
        quadLight_t light;
        
        light.a = loadVec3(arguments, 0);
        light.ab = loadVec3(arguments, 3);
        light.ac = loadVec3(arguments, 6);
        light.nl = cross(light.ab, light.ac);
        light.intensity = loadVec3(arguments, 9);
        
        _quadLights.push_back(light);
        
        // create 2 triangles for ray-light intersection

        material_t areaLightMaterial = {
            vec3(0.0f),         // diffuse
            vec3(0.0f),         // specular
            1.0f,               // shininess
            light.intensity,    // emission
            vec3(0.0f),         // ambient
            0.0f,               // roughness
            "phong"
        };

        vec3 a = light.a;
        vec3 b = light.a + light.ab;
        vec3 c = light.a + light.ac;
        vec3 d = light.a + light.ab + light.ac;

        loadTriangles(a, b, c, areaLightMaterial);
        loadTriangles(c, b, d, areaLightMaterial);

        // insert to table: (triangle id) -> (light id)
        int triID = _triMaterials.size() - 1;
        int lightID = _quadLights.size() - 1;
        _quadLightsTable[triID] = lightID;
        _quadLightsTable[triID - 1] = lightID;
     
    } else if (command == "lightsamples") {
        
        _lightsamples = stoi(arguments[0]);
        
    } else if (command == "lightstratify") {
        
        _lightstratify = (arguments[0] == "on") ? true : false;
        
    } else if (command == "attenuation") {

        _curAttenuation = loadVec3(arguments);

    } else if (command == "ambient") {

        _curMaterial.ambient = loadVec3(arguments);

    } else if (command == "diffuse") {

        _curMaterial.diffuse = loadVec3(arguments);

    } else if (command == "specular") {

        _curMaterial.specular = loadVec3(arguments);

    } else if (command == "shininess") {

        _curMaterial.shininess = stof(arguments[0]);

    } else if (command == "emission") {

        _curMaterial.emission = loadVec3(arguments);
        
    } else if (command == "roughness") {
        
        _curMaterial.roughness = stof(arguments[0]);
        
    } else if (command == "brdf") {
        
        _curMaterial.brdf = arguments[0];
        
    } else {

        cerr << "Unknown command in scene file: '" << command << "'" << endl;

    }
}

void SceneLoader::loadSceneData(const string& filePath)
{
    ifstream file(filePath);
    if (!file.is_open()) throw runtime_error("Could not open file: '" + filePath + "'");

    string line;
    while (getline(file, line)) {
        istringstream tokenStream(line);

        string command;
        tokenStream >> command;

        if (command.size() == 0 || command[0] == '#') continue;

        vector<string> arguments;
        string argument;
        while (tokenStream >> argument) {
            arguments.push_back(argument);
        }

        executeCommand(command, arguments);
    }
}

void SceneLoader::loadEmbreeTriangles(RTCScene embreeScene)
{
    RTCGeometry embreeTriangles = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);

    vec3* embreeVertices = reinterpret_cast<vec3*>(rtcSetNewGeometryBuffer(
        embreeTriangles,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        sizeof(vec3),
        _vertices.size()));
    memcpy(embreeVertices, _vertices.data(), _vertices.size() * sizeof(vec3));

    uvec3* embreeIndices = reinterpret_cast<uvec3*>(rtcSetNewGeometryBuffer(
        embreeTriangles,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        sizeof(uvec3),
        _indices.size()));
    memcpy(embreeIndices, _indices.data(), _indices.size() * sizeof(uvec3));

    rtcCommitGeometry(embreeTriangles);
    rtcAttachGeometry(embreeScene, embreeTriangles);
    rtcReleaseGeometry(embreeTriangles);
}

void SceneLoader::loadEmbreeSpheres(RTCScene embreeScene)
{
    RTCScene embreeSphereScene = rtcNewScene(_embreeDevice);

    RTCGeometry embreeSphere = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_SPHERE_POINT);

    vec4* embreeSpherePoint = reinterpret_cast<vec4*>(rtcSetNewGeometryBuffer(
        embreeSphere,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT4,
        sizeof(vec4),
        1));
    *embreeSpherePoint = vec4(0.0f, 0.0f, 0.0f, 1.0f);

    rtcCommitGeometry(embreeSphere);
    rtcAttachGeometry(embreeSphereScene, embreeSphere);
    rtcReleaseGeometry(embreeSphere);
    rtcCommitScene(embreeSphereScene);

    for (mat4 transform : _sphereTransforms) {
        RTCGeometry embreeSphereInstance = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_INSTANCE);
        rtcSetGeometryInstancedScene(embreeSphereInstance, embreeSphereScene);
        rtcSetGeometryTimeStepCount(embreeSphereInstance, 1);
        rtcSetGeometryTransform(
            embreeSphereInstance,
            0,
            RTC_FORMAT_FLOAT4X4_COLUMN_MAJOR,
            value_ptr(transform));
        rtcCommitGeometry(embreeSphereInstance);
        rtcAttachGeometry(embreeScene, embreeSphereInstance);
        rtcReleaseGeometry(embreeSphereInstance);
    }

    rtcReleaseScene(embreeSphereScene);
}

RTCScene SceneLoader::createEmbreeScene()
{
    RTCScene embreeScene = rtcNewScene(_embreeDevice);
    loadEmbreeTriangles(embreeScene);
    loadEmbreeSpheres(embreeScene);
    rtcCommitScene(embreeScene);
    return embreeScene;
}

Scene* SceneLoader::commitSceneData()
{
    float aspectRatio = static_cast<float>(_imageSize.x) / _imageSize.y;
    vec3 cameraLook = normalize(_cameraLookAt - _cameraOrigin);
    vec3 imagePlaneRight = normalize(cross(cameraLook, _cameraUp));
    vec3 imagePlaneUp = normalize(cross(imagePlaneRight, cameraLook));

    camera_t camera;
    camera.origin = _cameraOrigin;
    camera.imagePlaneTopLeft =
        _cameraOrigin
        + cameraLook / tan(PI * _cameraFieldOfView / 360.0f)
        + imagePlaneUp
        - aspectRatio * imagePlaneRight;

    camera.pixelRight = (2.0f * aspectRatio / _imageSize.x) * imagePlaneRight;
    camera.pixelDown = (-2.0f / _imageSize.y) * imagePlaneUp;

    vector<mat3> sphereNormalTransforms;
    for (size_t i = 0; i < _sphereTransforms.size(); i++) {
        sphereNormalTransforms.push_back(inverseTranspose(mat3(_sphereTransforms[i])));
    }

    Scene* scene = new Scene();
    scene->imageSize = _imageSize;
    scene->outputFileName = _outputFileName;
    scene->camera = camera;
    scene->sphereNormalTransforms = move(sphereNormalTransforms);
    scene->sphereMaterials = move(_sphereMaterials);
    scene->triMaterials = move(_triMaterials);
    scene->directionalLights = move(_directionalLights);
    scene->pointLights = move(_pointLights);
    scene->quadLights = move(_quadLights);
    scene->quadLightsTable = move(_quadLightsTable);
    scene->integrator = _integrator;
    scene->maxDepth = _maxDepth;
    scene->lightsamples = _lightsamples;
    scene->spp = _spp;
    scene->gamma = _gamma;
    scene->lightstratify = _lightstratify;
    scene->NEE = _NEE;
    scene->RR = _RR;
    scene->imptSampling = _imptSampling;
    scene->embreeScene = createEmbreeScene();

    return scene;
}

void loadScene(const string& filePath, RTCDevice embreeDevice, Scene** scene) 
{
    SceneLoader sceneLoader(embreeDevice);
    sceneLoader.loadSceneData(filePath);
    *scene = sceneLoader.commitSceneData();
}
