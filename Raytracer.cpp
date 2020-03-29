// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)

#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <string>
#include <map>
#include <list>
#define M_PI 3.1415926535897932

#include <random>
std::default_random_engine engine[8];
std::uniform_real_distribution<double> distrib(0, 1);

void save_image(const char* filename, const unsigned char* tableau, int w, int h)   // (0,0) is top-left corner
{

    FILE* f;

    int filesize = 54 + 3 * w * h;

    unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
    unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
    unsigned char bmppad[3] = { 0,0,0 };

    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[4] = (unsigned char)(w);
    bmpinfoheader[5] = (unsigned char)(w >> 8);
    bmpinfoheader[6] = (unsigned char)(w >> 16);
    bmpinfoheader[7] = (unsigned char)(w >> 24);
    bmpinfoheader[8] = (unsigned char)(h);
    bmpinfoheader[9] = (unsigned char)(h >> 8);
    bmpinfoheader[10] = (unsigned char)(h >> 16);
    bmpinfoheader[11] = (unsigned char)(h >> 24);

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    unsigned char* row = new unsigned char[w * 3];
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            row[j * 3] = tableau[(w * (h - i - 1) * 3) + j * 3 + 2];
            row[j * 3 + 1] = tableau[(w * (h - i - 1) * 3) + j * 3 + 1];
            row[j * 3 + 2] = tableau[(w * (h - i - 1) * 3) + j * 3];
        }
        fwrite(row, 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }
    fclose(f);
    delete[] row;
}

class Vector
{
public:
    Vector(double x = 0, double y = 0, double z = 0)
    {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    };
    double coord[3];

    double operator[](int i)  const
    {
        return coord[i];
    };

    double& operator[](int i)
    {
        return coord[i];
    };

    double getNorm2() const
    {
        return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
    }

    void Normalize()
    {
        double N = sqrt(getNorm2());
        coord[0] = coord[0] / N;
        coord[1] = coord[1] / N;
        coord[2] = coord[2] / N;
    }

    Vector getNormalized()
    {
        Vector result(*this);
        result.Normalize();
        return result;
    }

    Vector& operator+=(const Vector& A)
    {
        coord[0] += A[0];
        coord[1] += A[1];
        coord[2] += A[2];
        return *this;
    }
};

double dot(const Vector& A, const Vector& B)
{
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

Vector operator+(const Vector& A, const Vector& B)
{
    return Vector(A[0] + B[0], A[1] + B[1], A[2] + B[2]);
}

Vector operator-(const Vector& A, const Vector& B)
{
    return Vector(A[0] - B[0], A[1] - B[1], A[2] - B[2]);
}

Vector operator*(const double& a, const Vector& B)
{
    return Vector(a * B[0], a * B[1], a * B[2]);
}

Vector operator*(const Vector& A, const double& b)
{
    return Vector(A[0] * b, A[1] * b, A[2] * b);
}

Vector operator*(const Vector& A, const Vector& B)
{
    return Vector(A[0] * B[0], A[1] * B[1], A[2] * B[2]);
}

Vector operator/(const Vector& A, const double& b)
{
    return Vector(A[0] / b, A[1] / b, A[2] / b);
}

Vector cross(const Vector& A, const Vector& B)
{
    return Vector(A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]);
}

Vector random_cos(const Vector& N) {
    //orthogonal frame:
    Vector T1;
    if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) {
        T1 = Vector(0, -N[2], N[1]);
    }
    else if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])) {
        T1 = Vector(-N[2], 0, N[0]);
    }
    else {
        T1 = Vector(-N[1], N[0], 0);
    }
    T1.Normalize();
    Vector T2 = cross(N, T1);

    double r1 = distrib(engine[omp_get_thread_num()]);
    double r2 = distrib(engine[omp_get_thread_num()]);
    Vector V;
    V[0] = cos(2 * M_PI * r1) * sqrt(1 - r2);
    V[1] = sin(2 * M_PI * r1) * sqrt(1 - r2);
    V[2] = sqrt(r2);

    return V[0] * T1 + V[1] * T2 + V[2] * N;
}



class Ray
{
public:
    Ray() {};
    Ray(const Vector& C, Vector u) : C(C), u(u) {};
    Vector C, u;
};

class Object
{
public:
    Object() {};
    Object(const Vector& albedo, bool miroir = false, bool transparent = false, double emissivity = 0) : albedo(albedo), miroir(miroir), transparent(transparent), emissivity(emissivity) {};
    virtual bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& color) const = 0;

    Vector albedo;
    bool miroir;
    bool transparent;
    double emissivity;
};


class Sphere : public Object
{
public:
    Sphere(const Vector& O, double R, const Vector& albedo, bool miroir = false, bool transparent = false, double emissivity = 0) : O(O), R(R), Object(albedo, miroir, transparent, emissivity) {
    };

    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& color) const
    {
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = (r.C - O).getNorm2() - R * R;

        double delta = b * b - 4 * a * c;
        if (delta < 0)
            return false;
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);

        if (t2 < 0)
            return false;

        if (t1 > 0)
            t = t1;
        else
            t = t2;

        P = r.C + t * r.u;
        N = (P - O).getNormalized();
        color = albedo;
        return true;
    }

    Vector O;
    double R;

};


class Triangle : public Object 
{
public:
    Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& albedo, bool miroir = false, bool transparent = false, double emissivity = 0) : A(A), B(B), C(C), Object(albedo, miroir, transparent, emissivity) {
    };
    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& color) const {
        double beta, gamma;
        color = albedo;
        return intersect(r, P, N, t, beta, gamma);
    }

    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, double& beta, double& gamma) const {

        N = -1 * cross(B - A, C - A);
        N.Normalize();
        double num = -1 * dot(r.C - A, N);
        double den = dot(r.u, N);
        if (den == 0) return false;
        t = num / den;
        if (t < 0) {
            return false;
        }
        P = r.C + t * r.u;
        Vector CA = C - A;
        Vector BA = B - A;
        Vector PA = P - A;
        double dot00 = dot(CA,CA);
        double dot01 = dot(CA, BA);
        double dot02 = dot(CA, PA);
        double dot11 = dot(BA, BA);
        double dot12 = dot(BA, PA);
        double inv = 1 / (dot00 * dot11 - dot01 * dot01);
        gamma = (dot11 * dot02 - dot01 * dot12) * inv;
        beta = (dot00 * dot12 - dot01 * dot02) * inv;
        return(beta >= 0 && gamma >= 0 && (beta + gamma) < 1);
        

    }

    const Vector A, B, C;

};



class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
    };
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk;
    int ni, nj, nk;
    int faceGroup;
};


class BBox {
public:
    BBox() {};
    bool intersect(const Ray& r) const {

        bool rin = true;
        Vector pos, bmM;
        for (int i = 0; i < 3; i++) {
            if (r.C[i] < m[i]) {
                rin = false;
                pos[i] = -1;
                bmM[i] = m[i];
            }
            else if (r.C[i] > M[i]) {
                rin = false;
                pos[i] = 1;
                bmM[i] = M[i];
            }
            else {
                pos[i] = 0;
            }
        }
        if (rin) return true;

        Vector T;
        for (int i = 0; i < 3; i++) {
            if (pos[i] != 0 && r.u[i] != 0) {
                T[i] = (bmM[i] - r.C[i]) / r.u[i];
            }
            else {
                T[i] = -1;
            }
        }

        int planid = 0;
        for (int i = 1; i < 3; i++)
            if (T[planid] < T[i]) {
                planid = i;
            }
        if (T[planid] < 0) return (false);

        Vector P;
        for (int i = 0; i < 3; i++) {
            if (planid != i) {
                P[i] = r.C[i] + T[planid] * r.u[i];
                if (P[i] < m[i] || P[i] > M[i]) return (false);
            }
        }
        return true;
    }

    Vector m, M;
};


class BVH {
public:
    BVH* fg, * fd;
    int i0, i1;
    BBox bb;
};

class Geometry : public Object 
{
public:
    Geometry() {};
    Geometry(const char* obj, double scaling, const Vector& offset, const Vector& albedo, bool miroir = false, bool transparent = false, double emissivity = 0) : Object(albedo, miroir, transparent, emissivity) {

        readOBJ(obj);

        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * scaling + offset;
        }
        GetBVH(&bvh, 0, indices.size());


    }

    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");

        std::map<std::string, int> groupNames;
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                if (groupNames.find(std::string(grp)) != groupNames.end()) {
                    curGroup = groupNames[std::string(grp)];
                }
                else {
                    curGroup = groupNames.size();
                    groupNames[std::string(grp)] = curGroup;
                }
            }
            if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
                sscanf(line, "mtllib %[^\n]\n", matfile);
            }
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[2], &vec[1], &col[0], &col[1], &col[2]) == 6) {
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);  // helmet
                                                                                 //vec[2] = -vec[2]; //car2
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf_s(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;

                char* consumedline = line + 1;
                int offset;
                t.faceGroup = curGroup;
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;

                    indices.push_back(t);
                }
                else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }


                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.faceGroup = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }


        }
        fclose(f);


    }

    void add_texture(const char* filename) {

        textures.resize(textures.size() + 1);
        w.resize(w.size() + 1);
        h.resize(h.size() + 1);

        FILE* f;
        f = fopen(filename, "rb");
        unsigned char info[54];
        fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

        w[w.size() - 1] = *(int*)&info[18]; // extract image height and width from header
        h[h.size() - 1] = *(int*)&info[22];

        int size = 3 * w[w.size() - 1] * h[h.size() - 1];
        textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
        fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
        fclose(f);

        for (int i = 0; i < size; i += 3) {
            std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
        }
    }

    BBox GetBBox(int i0, int i1) {

        BBox b;
        b.m = vertices[indices[i0].vtxi];
        b.M = vertices[indices[i0].vtxi];
        for (int i = i0; i < i1; i++) {
            for (int j=0; j<3; j++) {
                b.m[j] = std::min(b.m[j], vertices[indices[i].vtxi][j]);
                b.M[j] = std::max(b.M[j], vertices[indices[i].vtxi][j]);
                b.m[j] = std::min(b.m[j], vertices[indices[i].vtxj][j]);
                b.M[j] = std::max(b.M[j], vertices[indices[i].vtxj][j]);
                b.m[j] = std::min(b.m[j], vertices[indices[i].vtxk][j]);
                b.M[j] = std::max(b.M[j], vertices[indices[i].vtxk][j]);
            }
        }
        return b;
    }

    void GetBVH(BVH* node, int i0, int i1) {
        node->bb = GetBBox(i0, i1);
        node->i0 = i0;
        node->i1 = i1;
        node->fg = NULL;
        node->fg = NULL;

        Vector diagonal = node->bb.M - node->bb.m;
        int dim;
        if ((diagonal[0] > diagonal[1]) && (diagonal[0] > diagonal[2])) {
            dim = 0;
        }
        else {
            if ((diagonal[1] > diagonal[0]) && (diagonal[1] > diagonal[2])) {
                dim = 1;
            }
            else {
                dim = 2;
            }
        }

        double split = node->bb.m[dim] + diagonal[dim] * 0.5;

        int pivot = i0;
        for (int i = i0; i < i1; i++) {
            double barycenter = (vertices[indices[i].vtxi][dim] + vertices[indices[i].vtxj][dim] + vertices[indices[i].vtxk][dim]) / 3;

            if (barycenter < split) {
                std::swap(indices[i], indices[pivot]);
                pivot++;
            }
        }

        if (i1 - i0 > 3 && pivot > i0 + 1 && pivot < i1) {
            node->fg = new BVH;
            node->fd = new BVH;
            GetBVH(node->fg, i0, pivot);
            GetBVH(node->fd, pivot, i1);
        }
    }

    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& color) const {
        
        if (!bvh.bb.intersect(r)) return false;
        bool has_inter = false;
        t = 1e99;

        std::list<const BVH*> l;
        l.push_front(&bvh);

        while (!l.empty()) {

            const BVH* current = l.front();
            l.pop_front();
            if (current->fg && current->fg->bb.intersect(r)) {
                l.push_back(current->fg);
            }
            if (current->fd && current->fd->bb.intersect(r)) {
                l.push_back(current->fd);
            }

            if (!current->fg) {

                for (int i = current->i0; i < current->i1; i++) {

                    Vector localP, localN;
                    double localt;
                    double beta, gamma;
                    bool tri_intersect = Triangle(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], albedo, miroir, transparent).intersect(r, localP, localN, localt, beta, gamma);
                    if (tri_intersect) {
                        has_inter = true;
                        if (localt < t) {
                            t = localt;
                            P = localP;
                            double alpha = 1 - beta - gamma;
                            N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
                            N.Normalize();
                            // texture
                            if (miroir == false && transparent == false && !textures.empty()) {
                                Vector UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                                int textureid = indices[i].faceGroup;
                                int Wi = w[textureid];
                                int Hi = h[textureid];
                                int x = UV[0] * (Wi - 1);
                                int y = UV[1] * (Hi - 1);

                                double cr = (textures[textureid][(y * Wi + x) * 3]) / 255.;
                                double cg = (textures[textureid][(y * Wi + x) * 3 + 1]) / 255.;
                                double cb = (textures[textureid][(y * Wi + x) * 3 + 2]) / 255.;
                                color = Vector(cr, cg, cb);
                            }
                            else {
                                color = Vector(1, 1, 1);
                            }
                            
                        }
                    }
                }

            }
        }

        return has_inter;
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vector> vertexcolors;

    std::vector<std::vector<unsigned char> > textures;
    std::vector<int> w, h;

    BVH bvh;
};


class Scene
{
public:
    Scene() {};
    void addSphere(Sphere& s)
    {
        objects.push_back(&s);
    }

    void addTriangle(Triangle& t)
    {
        objects.push_back(&t);
    }

    void addGeometry(Geometry& g)
    {
        objects.push_back(&g);
    }

    bool intersection(const Ray& r, Vector& P, Vector& N, int& sphere_id, Vector& color) const
    {

        bool has_inter = false;
        double min_t = 1E99;

        for (int i = 0; i < objects.size(); i++)
        {
            Vector localP, localN, localcolor;
            double t;
            bool local_has_inter = objects[i]->intersect(r, localP, localN, t, localcolor);
            if (local_has_inter)
            {
                has_inter = true;
                if (t < min_t)
                {
                    min_t = t;
                    P = localP;
                    N = localN;
                    color = localcolor;
                    sphere_id = i;
                }
            }
        }
        return has_inter;
    }

    std::vector<Object*> objects;
    Sphere *lumiere;
};


Vector getColor(const Ray& r, const Scene& s, int nb_rebond)
{

    Vector P, N, albedo;
    int sphere_id;
    bool has_inter = s.intersection(r, P, N, sphere_id, albedo);
    Vector intensite_pixel(0, 0, 0);

    if (has_inter)
    {
        if (nb_rebond == 0) return Vector(0, 0, 0);

        //miroir 
        if (s.objects[sphere_id]->miroir)
        {
            Vector direction_refl = r.u - 2 * dot(r.u, N) * N;
            Ray reflechi(P + 0.001 * N, direction_refl);
            intensite_pixel = getColor(reflechi, s, nb_rebond - 1);
        }
        else
            // transparence
            if (s.objects[sphere_id]->transparent)
            {
                double n1 = 1;
                double n2 = 1.45;
                Vector N_trans(N);
                Ray new_r;
                if (dot(r.u, N) > 0) {
                    n1 = 1.45;
                    n2 = 1;
                    N_trans = -1 * N;
                }
                double radical = 1 - pow(n1 / n2, 2) * (1 - pow(dot(N_trans, r.u), 2));
                if (radical > 0) {
                    Vector direction_refr = (n1 / n2) * (r.u - dot(r.u, N_trans) * N_trans) - N_trans * sqrt(radical);
                    direction_refr.Normalize();

                    double R0 = pow((n1 - n2) / (n1 + n2), 2);
                    double R;
                    if (n1<n2) {
                        R = R0 + (1 - R0) * pow(1 + dot(r.u, N), 5);
                    }
                    else {
                        R = R0 + (1 - R0) * pow(1 - dot(direction_refr, N), 5);
                    }
                    double r1 = distrib(engine[omp_get_thread_num()]);

                    if (r1 < R) {
                        Vector direction_refl = r.u - 2 * dot(r.u, N) * N;
                        direction_refl.Normalize();
                        new_r = Ray(P + 0.001 * N_trans, direction_refl);
                    }
                    else {
                        new_r = Ray(P - 0.001 * N_trans, direction_refr);
                    }
                }
                else {
                    Vector direction_refl = r.u - 2 * dot(r.u, N_trans) * N_trans;
                    direction_refl.Normalize();
                    new_r = Ray(P + 0.001 * N_trans, direction_refl);
                }
                intensite_pixel = getColor(new_r, s, nb_rebond - 1);
            }

            else
            {
                // Contribution directe
                Vector axeOP = (P - s.lumiere->O).getNormalized();
                Vector dir_aleatoire = random_cos(axeOP);
                dir_aleatoire.Normalize();
                Vector point_aleatoire = dir_aleatoire * s.lumiere->R + s.lumiere->O;
                Vector wi = (point_aleatoire - P).getNormalized();
                double d_light = (point_aleatoire - P).getNorm2();
                Vector Np = dir_aleatoire;

                double costheta = std::max(0., dot(N, wi));
                double costhetaprime = dot(dir_aleatoire, -1 * wi);
                double costhetaseconde = dot(axeOP, dir_aleatoire);

                Ray ray_light(P + 0.01 * N, wi);
                Vector P_light, N_light, albedo_light;
                int sphere_id_light;
                bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, albedo_light);
                if (has_inter_light && (P_light - P).getNorm2() < d_light-0.01) {
                    intensite_pixel = Vector(0, 0, 0);
                }
                else {
                    intensite_pixel = (s.lumiere->emissivity / (4 * M_PI * d_light) * costheta * costhetaprime / costhetaseconde) * albedo;
                }
                

                // Contribution indirecte
                Vector direction_aleatoire = random_cos(N);
                Ray rayon_aleatoire(P + 0.001 * N, direction_aleatoire);

                intensite_pixel += getColor(rayon_aleatoire, s, nb_rebond - 1) * albedo;

            }
    }
    return intensite_pixel;
}

int main()
{
    int W = 1024;
    int H = 1024;
    int Nrays = 100;
    double fov = 60 * M_PI / 180;

    // source de lumière
    Sphere slum(Vector(-15, 20, 40), 2, Vector(1., 1., 1.), false, false, 6e9);

    // sphère principale
    Sphere sphere1(Vector(0, 0, 0), 10, Vector(1, 1, 1),false,false);

    // test de transparence
    Sphere sphere2(Vector(20, 5, 0), 10, Vector(1, 0, 0), false, true);

    // test sphère miroir
    Sphere sphere3(Vector(-20, 5, 0), 10, Vector(1, 0, 0), true, false);

    // perspective

    Sphere sphere4(Vector(10, 5, 40), 5, Vector(1, 1, 1));

    // décor
    Sphere s2(Vector(0, -1000, 0), 990, Vector(0.043, 0.043, 0.043)); // sol
    Sphere s3(Vector(0, 1000, 0), 940, Vector(1, 1, 1)); // plafond
    Sphere s4(Vector(-990, 0, 0), 940, Vector(0, 1, 0)); // mur droit
    Sphere s5(Vector(990, 0, 0), 940, Vector(0, 0, 1)); // mur gauche
    Sphere s6(Vector(0, 0, -1000), 940, Vector(0, 1, 1)); // mur fond

    // Triangle

    Triangle tri(Vector(20,0,0), Vector(-20,0,0), Vector(0,30,0),Vector(1,0,0));

    // Maillage
    Geometry g1 =  Geometry("Girl.obj", 15, Vector(0, -10, 20), Vector(1, 1, 1),false,false);
    g1.add_texture("12c14c70.bmp");
    g1.add_texture("13932ef0.bmp");
    g1.add_texture("19d89130.bmp");
    g1.add_texture("16cecd10.bmp");
    g1.add_texture("16c2e0d0.bmp");
    g1.add_texture("12dbd6d0.bmp");


    // Création de scène
    Scene s;
    
    s.addSphere(slum);
 
    s.addSphere(sphere1);
    //s.addSphere(sphere2);
    //s.addSphere(sphere3);
    //s.addSphere(sphere4);
    //s.addGeometry(g1);
    //s.addTriangle(tri);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);
    s.lumiere = &slum;

    Vector position_camera(0, 0, 55);
    double focus_distance = 55;


    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {

            Vector color(0, 0, 0);
            for (int k = 0; k < Nrays; k++) {

                //Anti-aliasing
                double r1 = distrib(engine[omp_get_thread_num()]);
                double r2 = distrib(engine[omp_get_thread_num()]);
                double R = sqrt(-2 * log(r1));
                double dx = R * cos(2 * M_PI * r2);
                double dy = R * sin(2 * M_PI * r2);

                Vector u(-j + (W / 2) + 0.5 + dx, -i + (H / 2) + 0.5 + dy, -W / (2 * tan(fov / 2)));
                u.Normalize();

                // Profondeur de champs
                double rp = distrib(engine[omp_get_thread_num()]);
                double dx_aperture = R * cos(2 * M_PI * rp) * 0.25;
                double dy_aperture = R * sin(2 * M_PI * rp) * 0.25;

                Vector destination = position_camera + focus_distance * u;
                Vector new_origin = position_camera + Vector(dx_aperture, dy_aperture, 0);

                Ray r(new_origin, (destination - new_origin).getNormalized());

                color += getColor(r, s, 5) / Nrays;
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 0.45)));
            image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 0.45)));
            image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 0.45)));

        }
    }

    //    img[(100 * W + 500)*3] = 255;  // pixel at (x, y) = (50, 10) is red

    save_image("test.bmp", &image[0], W, H);

    return 0;
}
