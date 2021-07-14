#pragma once

#include <glm/glm.hpp>

#include "Scene.h"
#include "Constants.h"

using namespace glm;
using namespace std;

class Integrator {

protected:

    Scene* _scene;
    
public:

    void setScene(Scene* scene) { _scene = scene; }

    virtual vec3 traceRay(vec3& origin, vec3& direction, vec3& target) = 0;

};


class RayTracerIntegrator : public Integrator {

private:

    float getRand();

    float dotProd(vec3& a, vec3& b);

    vec3 computeShading(const vec3& incidentDirection, const vec3& toLight, const vec3& normal,
                        const vec3& brightness, const material_t& material);
    
    vec3 getIrradVector(vector<vec3>& v, vec3& r);
    
    vector<vec3> getLightSamplePoints(quadLight_t& light, float N);

    vec3 rotateToCoord(vec3& s, vec3& z);
        
    vec3 getHemisphereSample(vec3& n);
    
    vec3 getCosineSample(vec3& n);
    
    vec3 getBRDFSample(vec3& z, float shininess);
    
    vec3 getBRDFHalfVector(vec3& n, float roughness);
        
    float microfacet(float a, float hn);
    
    float Smith(float a, float vn);
    
    vec3 Schlick(vec3& ks, float wih);
    
    bool rayIntersectLight(vec3& p0, vec3& p1, vec3& xp, quadLight_t& light);
    
    bool rayIntersectTri(vec3& A, vec3& B, vec3& C, vec3& p0, vec3& p1, vec3& xp);
        
    float pdfNEE(vec3& x, vec3& wi);

    vec3 getSampleDirection(vec3& n, vec3& direction, material_t& m, float& t);
        
    vec3 getBRDF(vec3& n, vec3& direction, vec3& wi, material_t& m);
    
    float getPDF(vec3& n, vec3& direction, vec3& wi, material_t& m, float t);
    
    vec3 getThroughput(vec3& n, vec3& direction, vec3& wi, material_t& m, float t);

    vec3 rayTracing(vec3& origin, vec3& direction, int depth);

    vec3 analyticDirectLighting(vec3& origin, vec3& direction);

    vec3 directLighting(vec3& origin, vec3& direction, bool emission);

    vec3 directLightingMIS(vec3& origin, vec3& direction);
        
    vec3 pathTracingBasic(vec3& origin, vec3& direction, int depth);
    
    vec3 pathTracingNEE(vec3& origin, vec3& direction, vec3& T, int depth);
        
    vec3 pathTracing(vec3& origin, vec3& target);

public:

    virtual vec3 traceRay(vec3& origin, vec3& direction, vec3& target);

};


