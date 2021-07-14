#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtc/random.hpp>

#include "Integrator.h"

float RayTracerIntegrator::getRand() {
    return (float) rand() / RAND_MAX;
}

float RayTracerIntegrator::dotProd(vec3& a, vec3& b) {
    return clamp(dot(a, b), 0.0f, 1.0f);
}

vec3 RayTracerIntegrator::computeShading(const vec3& incidentDirection, const vec3& toLight, const vec3& normal,
                                         const vec3& brightness, const material_t& material) {
    
    vec3 diffuseReflectance = material.diffuse * fmax(dot(normal, toLight), 0.0f);
    vec3 halfAngle = normalize(toLight + incidentDirection);
    vec3 specularReflectance = material.specular * pow(fmax(dot(normal, halfAngle), 0.0f), material.shininess);
    return brightness * (diffuseReflectance + specularReflectance);
    
}

vec3 RayTracerIntegrator::rayTracing(vec3& origin, vec3& direction, int depth) {

    vec3 color = vec3(0.0f);

    vec3 x;         // hit posision
    vec3 n;         // hit normal
    material_t m;   // hit material
    
    bool hit = _scene->castRay(origin, direction, &x, &n, &m);
    vec3 oppositeDir = -direction;
    
    if (hit) {

        color += m.ambient + m.emission;

        for (const directionalLight_t light : _scene->directionalLights) {
            bool occluded = _scene->castOcclusionRay(x, light.toLight);
            if (!occluded) {
                color += computeShading(oppositeDir, light.toLight, n, light.brightness, m);
            }
        }
        for (const pointLight_t light : _scene->pointLights) {
            vec3 toLight = light.point - x;
            float lightDistance = length(toLight);
            toLight /= lightDistance;

            bool occluded = _scene->castOcclusionRay(x, toLight, lightDistance);
            if (!occluded) {
                float falloff = light.attenuation.x +
                                light.attenuation.y * lightDistance +
                                light.attenuation.z * lightDistance * lightDistance;
                color += computeShading(oppositeDir, toLight, n, light.brightness / falloff, m);
            }
        }

        // trace reflected ray if not reaching max depth 
        if (depth < _scene->maxDepth) {
            vec3 r = reflect(direction, n);
            color += m.specular * rayTracing(x, r, depth + 1);
        }
    }
    return color;
}

vec3 RayTracerIntegrator::getIrradVector(vector<vec3>& v, vec3& r) {
    
    vec3 irradVec = vec3(0.0f);
    unsigned int n = v.size();
    
    for (unsigned int k = 0; k < n; k++) {
        vec3 a = v[k % n] - r, b = v[(k + 1) % n] - r;
        float theta = acos(dot(a / length(a), b / length(b)));
        
        vec3 gamma = cross(a, b);
        gamma /= length(gamma);
        
        irradVec += theta * gamma;
    }
    
    return 0.5f * irradVec;
    
}

vec3 RayTracerIntegrator::analyticDirectLighting(vec3& origin, vec3& direction) {
    
    vec3 outputColor = vec3(0.0f);

    vec3 x;         // hit posision
    vec3 n;         // hit normal
    material_t m;   // hit material
    
    bool hit = _scene->castRay(origin, direction, &x, &n, &m);
    if (!hit) return outputColor;

    for (const quadLight_t light : _scene->quadLights) {
        vector<vec3> v = {light.a, light.a + light.ab, light.a + light.ab + light.ac, light.a + light.ac};
        vec3 irradVec = getIrradVector(v, x);
        outputColor += m.emission + m.diffuse * INV_PI * light.intensity * dot(irradVec, n);
    }
    
    return outputColor;
}

vector<vec3> RayTracerIntegrator::getLightSamplePoints(quadLight_t& light, float N) {
    
    bool stratify = _scene->lightstratify;
    vector<vec3> samplePoints;
    
    if (!stratify) { // random sampling points
        for (int k = 0; k < N; k++) {
            float u1 = getRand(), u2 = getRand();
            vec3 point = light.a + u1 * light.ab + u2 * light.ac;
            samplePoints.push_back(point);
        }
    }
    else { // stratified sampling points
        int M = (int) sqrt(N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                float u1 = getRand(), u2 = getRand();
                vec3 point = light.a + ((j + u1) / M) * light.ab + ((i + u2) / M) * light.ac;
                samplePoints.push_back(point);
            }
        }
    }

    return samplePoints;
}

vec3 RayTracerIntegrator::directLighting(vec3& origin, vec3& direction, bool emission) {
    
    vec3 color = vec3(0.0f);

    vec3 x;         // hit posision
    vec3 n;         // hit normal
    material_t m;   // hit material
    
    bool hit = _scene->castRay(origin, direction, &x, &n, &m);
    if (!hit) return color;
    
    int N = _scene->lightsamples;
    
    for (quadLight_t light : _scene->quadLights) {
        
        vec3 nl = light.nl;
        float A = length(nl);
        nl = nl / A;
        
        // generate sample points on area light
        vector<vec3> samplePoints = getLightSamplePoints(light, N);
        
        vec3 sum(0.0f);
        for (auto xp : samplePoints) {
            
            vec3 wi = normalize(xp - x);
            float R = length(xp - x);
            
            bool occluded = _scene->castOcclusionRay(x, wi, R);
            if (occluded) continue;
                
            vec3 f = getBRDF(n, direction, wi, m);
            float G = fmax(0.0f, dot(n, wi)) * fmax(0.0f, dot(nl, wi)) / (R * R);
            sum += f * G;
        }
        
        vec3 Le = emission ? m.emission : vec3(0.0f);
        color += Le + light.intensity * (A / (float) N) * sum;
    }
    
    return color;
}

vec3 RayTracerIntegrator::rotateToCoord(vec3& s, vec3& z) {

    vec3 w = z;
    vec3 a = sphericalRand(1.0f);
    vec3 u = normalize(cross(a, w));
    vec3 v = cross(w, u);
    
    return s.x * u + s.y * v + s.z * w;
}

vec3 RayTracerIntegrator::getHemisphereSample(vec3& n) {
        
    float u1 = getRand();
    float u2 = getRand();
    
    float a1 = acos(u1);
    float a2 = TWO_PI * u2;
    
    vec3 s = vec3(cos(a2) * sin(a1),
                  sin(a2) * sin(a1),
                  cos(a1));
    
    return rotateToCoord(s, n);
}

vec3 RayTracerIntegrator::getCosineSample(vec3& n) {
    
    float u1 = getRand();
    float u2 = getRand();
    
    float a1 = acos(sqrt(u1));
    float a2 = TWO_PI * u2;
    
    vec3 s = vec3(cos(a2) * sin(a1),
                  sin(a2) * sin(a1),
                  cos(a1));
    
    return rotateToCoord(s, n);
}

vec3 RayTracerIntegrator::getBRDFSample(vec3& z, float shininess) {
    
    float u1 = getRand();
    float u2 = getRand();
    
    float a1 = acos(powf(u1, 1.0f / (shininess + 1.0f)));
    float a2 = TWO_PI * u2;
    
    vec3 s = vec3(cos(a2) * sin(a1),
                  sin(a2) * sin(a1),
                  cos(a1));
    
    return rotateToCoord(s, z);
}

vec3 RayTracerIntegrator::getBRDFHalfVector(vec3& n, float roughness) {
    
    float u1 = getRand();
    float u2 = getRand();
    
    float a1 = atan(roughness * sqrt(u1) / sqrt(1.0f - u1));
    float a2 = TWO_PI * u2;
    
    vec3 s = vec3(cos(a2) * sin(a1),
                  sin(a2) * sin(a1),
                  cos(a1));
    
    return rotateToCoord(s, n);

}

vec3 RayTracerIntegrator::getSampleDirection(vec3& n, vec3& direction, material_t& m, float& t) {
    
    if (_scene->imptSampling == "cosine") {
        return getCosineSample(n);

    } else if (_scene->imptSampling == "brdf") {
        
        vec3 kd = m.diffuse;
        vec3 ks = m.specular;
        float s = m.shininess;
        float a = m.roughness;
        float kdGray = (kd.x + kd.y + kd.z) / 3.0f;
        float ksGray = (ks.x + ks.y + ks.z) / 3.0f;
        
        if (m.brdf == "ggx") { // ggx sample
            
            t = 1.0f; // when kd = ks = 0, sample Frensel reflection (t = 1)
            if (kdGray + ksGray != 0.0f) {
                t = fmax(0.25f, ksGray / (kdGray + ksGray));
            }
                        
            if (getRand() < t) { // specular
                vec3 h = getBRDFHalfVector(n, a);
                return reflect(direction, h);
            } else { // diffuse
                return getBRDFSample(n, 1);
            }
        } else { // phong sample
            
            t = ksGray / (kdGray + ksGray);
            
            if (getRand() < t) { // specular
                vec3 r = reflect(direction, n);
                return getBRDFSample(r, s);
            } else { // diffuse
                return getBRDFSample(n, 1);
            }
        }
    } else {
        return getHemisphereSample(n);
    }
}

vec3 RayTracerIntegrator::getBRDF(vec3& n, vec3& direction, vec3& wi, material_t& m) {
    
    vec3 kd = m.diffuse;
    vec3 ks = m.specular;
    float a = m.roughness;
    float s = m.shininess;
    
    if (m.brdf == "ggx") { // ggx
        
        vec3 w0 = -direction;
        
        vec3 fggx = vec3(0.0f);
        
        if (dot(w0, n) > 0.0f && dot(wi, n) > 0.0f) {

            vec3 h = normalize(w0 + wi);
    
            float win = dotProd(wi, n);
            float w0n = dotProd(w0, n);
            float wih = dotProd(wi, h);
            float hn = dotProd(h, n);
            
            vec3 F = Schlick(ks, wih);
            float G = Smith(a, win) * Smith(a, w0n);
            float D = microfacet(a, hn);
            fggx = F * G * D / (4.0f * win * w0n);
        }
        return kd * INV_PI + fggx;

    } else { // phong

        vec3 r = reflect(direction, n);
        return kd * INV_PI + ks * powf(dot(r, wi), s) * (s + 2.0f) / TWO_PI;
    }
}


float RayTracerIntegrator::getPDF(vec3& n, vec3& direction, vec3& wi, material_t& m, float t) {
    
    float pdf = 0.0f;
    float win = dot(n, wi);
    
    if (_scene->imptSampling == "cosine") { // cosine sampling
        
        pdf = win * INV_PI;
        
    } else if (_scene->imptSampling == "brdf") { // brdf sampling

        float s = m.shininess;
        float a = m.roughness;

        if (m.brdf == "ggx") { // ggx brdf

            vec3 w0 = -direction;
            vec3 h = normalize(w0 + wi);
            float hn = dotProd(h, n);
            pdf = (1.0f - t) * win * INV_PI +
                  t * (microfacet(a, hn) * dot(n, h) / (4.0f * dot(wi, h)));

        } else { // phong brdf

            vec3 r = reflect(direction, n);
            pdf = (1.0f - t) * win * INV_PI + t * pow(dot(r, wi), s) * (s + 1.0f) * INV_TWO_PI;
        }
        
    } else { // hemisphere sampling
        
        pdf = 1.0f / TWO_PI;
    }
    return pdf;
}

vec3 RayTracerIntegrator::getThroughput(vec3& n, vec3& direction, vec3& wi, material_t& m, float t) {
    
    float win = dotProd(wi, n);
    vec3 f = getBRDF(n, direction, wi, m);
    float pdf = getPDF(n, direction, wi, m, t);
    return f * win / pdf;
}

float RayTracerIntegrator::microfacet(float a, float hn) {
    
    float th = acos(hn);
    return a * a / (PI * powf(cos(th), 4.0f) * powf(a * a + powf(tan(th), 2.0f), 2.0f));
}

float RayTracerIntegrator::Smith(float a, float vn) {
    
    float tv = acos(vn);
    return 2.0f / (1.0f + sqrt(1.0f + a * a * powf(tan(tv), 2.0f)));
}

vec3 RayTracerIntegrator::Schlick(vec3& ks, float wih) {
    
    return ks + powf(1.0f - wih, 5.0f) * (vec3(1.0f) - ks);
}

bool RayTracerIntegrator::rayIntersectLight(vec3& p0, vec3& p1, vec3& xp, quadLight_t& light) {

    vec3 a = light.a;
    vec3 b = light.a + light.ab;
    vec3 c = light.a + light.ac;
    vec3 d = light.a + light.ab + light.ac;

    if (dot(p1, light.nl) < 0.0f) return false;
    
    if (rayIntersectTri(a, b, c, p0, p1, xp) ||
        rayIntersectTri(c, b, d, p0, p1, xp)) {
        return true;
    }

    return false;
}

bool RayTracerIntegrator::rayIntersectTri(vec3& A, vec3& B, vec3& C, vec3& p0, vec3& p1, vec3& xp) {

    float epsilon = 0.0001f;

    vec3 e1 = B - A;
    vec3 e2 = C - A;
    vec3 h, s, q;
    float a, f, u, v;

    h = cross(p1, e2);
    a = dot(e1, h);
    if (a > -epsilon && a < epsilon) return false;

    f = 1.0f / a;
    s = p0 - A;
    u = f * dot(s, h);
    if (u < 0.0f || u > 1.0f) return false;

    q = cross(s, e1);
    v = f * dot(p1, q);
    if (v < 0.0f || u + v > 1.0f) return false;

    float t = f * dot(e2, q);
    if (t < epsilon) return false;

    xp = p0 + t * p1;
    return true;
}

float RayTracerIntegrator::pdfNEE(vec3& x, vec3& wi) {
    
    float NEEpdf = 0.0f;
    int N = _scene->quadLights.size();
    
    for (quadLight_t light : _scene->quadLights) {
        
        vec3 xp; // hit position on light
        if (rayIntersectLight(x, wi, xp, light)) {
            float R = length(xp - x);
            vec3 nl = light.nl;
            float A = length(nl);
            nl = nl / A;
            NEEpdf += (R * R) / (A * abs(dot(nl, wi)));
        }
    }
    
    return (1.0f / (float) N) * NEEpdf;
}

vec3 RayTracerIntegrator::pathTracingBasic(vec3& origin, vec3& direction, int depth) {
    
    vec3 color = vec3(0.0f);
    
    vec3 x;         // hit posision
    vec3 n;         // hit normal
    material_t m;   // hit material
    
    bool hit = _scene->castRay(origin, direction, &x, &n, &m);
    if (!hit) return color;
        
    if (depth == 0 && length(m.emission) > 0.0f) return m.emission;
    if (depth == _scene->maxDepth) return color;
            
    float t = 0.0f;
    
    vec3 wi = getSampleDirection(n, direction, m, t);
    vec3 L = pathTracingBasic(x, wi, depth + 1);
    vec3 T = getThroughput(n, direction, wi, m, t);

    color = m.emission + T * L;
        
    return color;
}

vec3 RayTracerIntegrator::pathTracingNEE(vec3& origin, vec3& direction, vec3& T, int depth) {
    
    vec3 color = vec3(0.0f);
    
    if (depth == 0) {
        
        if (_scene->NEE == "on") { // NEE
            color += directLighting(origin, direction, true);
        } else { // MIS
            color += directLightingMIS(origin, direction);
        }
        color += pathTracingNEE(origin, direction, T, 1);
        return color;
    }
    
    if (depth == _scene->maxDepth) return color;

    // russian roulette
    float q = 0.0f;
    if (_scene->RR) {
        q = 1.0f - fmin(fmax(fmax(T.x, T.y), T.z), 1.0f);
        if (getRand() < q) return color;
    }
    
    vec3 x;         // hit posision
    vec3 n;         // hit normal
    material_t m;   // hit material
    
    bool hit = _scene->castRay(origin, direction, &x, &n, &m);
    if (!hit || length(m.emission) > 0.0f) return color;
    
    float t = 0.0f;
    vec3 wi = getSampleDirection(n, direction, m, t);
    vec3 newT = getThroughput(n, direction, wi, m, t) / (1.0f - q);
    vec3 Ld = directLighting(x, wi, false);
    vec3 Li = pathTracingNEE(x, wi, newT, depth + 1);

    color = newT * (Ld + Li);
    
    return color;
}

vec3 RayTracerIntegrator::directLightingMIS(vec3& origin, vec3& direction) {
    
    vec3 color = vec3(0.0f);

    vec3 x, xp;         // hit posision
    vec3 n, np;         // hit normal
    material_t m, mp;   // hit material
    
    bool hit = _scene->castRay(origin, direction, &x, &n, &m);
    if (!hit || length(m.emission) > 0.0f) return color;
    
    float t = 0.0f;
    
    // get the brdf sample direction wi
    vec3 wi = getSampleDirection(n, direction, m, t);
    
    // find if sample direction hit light source
    int hitID;
    hit = _scene->castRay(x, wi, &xp, &np, &mp, &hitID);

    if (hit && _scene->quadLightsTable.count(hitID) > 0) {

        int lightID = _scene->quadLightsTable[hitID];
        quadLight_t light = _scene->quadLights[lightID];
        if (dot(wi, light.nl) > 0.0f) {
            float NEEpdf = pdfNEE(x, wi);
            float BRDFpdf = getPDF(n, direction, wi, m, t);
            float weight = powf(BRDFpdf, 2.0f) / (powf(BRDFpdf, 2.0f) + powf(NEEpdf, 2.0f));
            vec3 f = getBRDF(n, direction, wi, m);
            float win = dotProd(wi, n);
            vec3 L = light.intensity;
            
            color += weight * L * f * win / BRDFpdf;
        }
    }
    
    float factor = (1.0f / (float) _scene->quadLights.size());

    for (quadLight_t light : _scene->quadLights) {

        // generate one sample point on area light
        vector<vec3> samplePoints = getLightSamplePoints(light, 1);
        vec3 xp = samplePoints[0];
        wi = normalize(xp - x);

        // test visibility
        bool occluded = _scene->castOcclusionRay(x, wi, length(xp - x));
        if (occluded) continue;

        float NEEpdf = pdfNEE(x, wi);
        if (NEEpdf == 0.0f) continue;
        float BRDFpdf = getPDF(n, direction, wi, m, t);
        float weight = powf(NEEpdf, 2.0f) / (powf(BRDFpdf, 2.0f) + powf(NEEpdf, 2.0f));
        vec3 f = getBRDF(n, direction, wi, m);
        float win = dotProd(wi, n);
        vec3 L = light.intensity;
        
        color += factor * weight * L * f * win / NEEpdf;
    }
    
    return color;
}

vec3 RayTracerIntegrator::pathTracing(vec3& origin, vec3& target) {
    
    vec3 color = vec3(0.0f);
    float weight = 1.0f / (float) _scene->spp;
    
    for (int i = 0; i < _scene->spp; i++) {
                    
        // u1, u2 rand in [-0.5, 0.5]
        float u1 = getRand() - 0.5f;
        float u2 = getRand() - 0.5f;
        
        // place the first sample in the center of the pixel
        if (i == 0) u1 = u2 = 0.0f;
        
        vec3 sampleTarget = target + u1 * _scene->camera.pixelRight + u2 * _scene->camera.pixelDown;
        vec3 sampleDir = normalize(sampleTarget - origin);
        
        if (_scene->NEE == "off") { // simple path tracer
            color += weight * pathTracingBasic(origin, sampleDir, 0);
        } else { // NEE
            vec3 T = vec3(1.0f);
            color += weight * pathTracingNEE(origin, sampleDir, T, 0);
        }
    }
    return color;
}

vec3 RayTracerIntegrator::traceRay(vec3& origin, vec3& direction, vec3& target) {
    
    vec3 color = vec3(0.0f);
    
    if (_scene->integrator == "analyticdirect") {
        color = analyticDirectLighting(origin, direction);
    }
    else if (_scene->integrator == "direct") {
        color = directLighting(origin, direction, true);
    }
    else if (_scene->integrator == "pathtracer") {
        color = pathTracing(origin, target);
    }
    else {
        color = rayTracing(origin, direction, 0);
    }
    
    color.x = powf(color.x, 1.0f / _scene->gamma);
    color.y = powf(color.y, 1.0f / _scene->gamma);
    color.z = powf(color.z, 1.0f / _scene->gamma);
    
    return color;
}
