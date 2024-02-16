#include <bits/stdc++.h>
#include <GL/glut.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "1805100_Point.hpp"
#include "1805100_Vec.hpp"
#include "1805100_BitmapImage.hpp"

#define _USE_MATH_DEFINES

using namespace std;

// Type aliases
using Point = Point3D<double>;
using Vec3 = Vec<double, 3>;

// Debug utility
#ifdef DEBUG
#define dbg(args...) { debugPrint(#args, args); }
template <typename... Args>
void debugPrint(const std::string& arg_names, const Args&... args) {
    stringstream ss;
    ss << ">> " << arg_names << " = ";
    ((ss << args << ", "), ...);
    string s(ss.str());
    s.erase(s.end()-2);
    std::cerr << s << std::endl;
}
#else
#define dbg(args...) (void) "Why trace ray?"
#endif

// Constants
const double EPSILON = 1e-10;
const double ROTATION_ANGLE_DELTA = 10.0;
const double CAMERA_DELTA = 10.0;
const double TILT_ANGLE_DELTA = 5.0;


// Utility functions
inline double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

inline double RadToDeg(double rad) {
    return rad * 180.0 / M_PI;
}

inline bool isEqual(double a, double b) {
    return abs(a - b) <= EPSILON;
}

struct Camera {
    Point pos;
    Vec3 look, right, up;

    Camera(const Point& pos, const Vec3& look, const Vec3& right, const Vec3& up) :
        pos(pos), look(look), right(right), up(up) 
    {
        this->look.normalize();
        this->right.normalize();
        this->up.normalize();
    }

    // Rotate the camera left around its up vector
    void rotateLeft(float angleDegrees) {
        Vec3 newLook = look.rotate(up, angleDegrees);
        Vec3 newRight = right.rotate(up, angleDegrees);
        look = newLook;
        right = newRight;
    }

    // Rotate the camera right around its up vector
    void rotateRight(float angleDegrees) {
        rotateLeft(-angleDegrees);
    }

    // Rotate the camera up around its right vector
    void rotateUp(float angleDegrees) {
        Vec3 newLook = look.rotate(right, angleDegrees);
        Vec3 newUp = up.rotate(right, angleDegrees);
        look = newLook;
        up = newUp;
    }

    // Rotate the camera down around its right vector
    void rotateDown(float angleDegrees) {
        rotateUp(-angleDegrees);
    }

    // Tilt the camera clockwise around its look vector
    void tiltClockwise(float angleDegrees) {
        Vec3 newRight = right.rotate(look, angleDegrees);
        Vec3 newUp = up.rotate(look, angleDegrees);
        right = newRight;
        up = newUp;
    }

    // Tilt the camera counterclockwise around its right vector
    void tiltCounterclockwise(float angleDegrees) {
        tiltClockwise(-angleDegrees);
    }

    // Move the camera left along its right vector
    void moveLeft(double distance) {
        pos = pos - right * distance;
    }

    // Move the camera right along its right vector
    void moveRight(double distance) {
        moveLeft(-distance);
    }

    // Move the camera up along its up vector
    void moveUp(double distance) {
        pos = pos + up * distance;
    }

    // Move the camera down along its up vector
    void moveDown(double distance) {
        moveUp(-distance);
    }

    // Move the camera forward along its look vector
    void moveForward(double distance) {
        pos = pos + look * distance;
    }

    // Move the camera backward along its look vector
    void moveBackward(double distance) {
        moveForward(-distance);
    }

    // Set up the camera using gluLookAt
    void setupCamera() const {
        gluLookAt(
            pos.x, pos.y, pos.z,
            pos.x + look[0], pos.y + look[1], pos.z + look[2],
            up[0], up[1], up[2]
        );
    }
};

Camera camera(
    Point(0, -100, 100),
    Vec3({0, 1, 0}),
    Vec3({1, 0, 0}),
    Vec3({0, 0, 1})
);

struct Ray {
    Point src;
    Vec3 dir;

    Ray(const Point& src, const Vec3& dir) : 
        src(src), dir(dir)
    {
        this->dir.normalize();
    }

    Ray getReflection(const Point& point, const Vec3& normal) const {
        Vec3 reflectedDir = dir - normal * (2 * dir.dot(normal));
        reflectedDir.normalize();
        return Ray(point, reflectedDir);
    }

    friend ostream& operator<<(ostream& out, const Ray& ray) {
        out << ray.src << " -> " << ray.dir;
        return out;
    }
};

struct ColorProp {
    double r, g, b;
    double ambient, diffuse, specular, reflection, shininess;

    ColorProp() :
        r(0), g(0), b(0),
        ambient(0),
        diffuse(0),
        specular(0),
        reflection(0),
        shininess(0)
    { }

    ColorProp(
        double r, double g, double b, 
        double ambient, 
        double diffuse, 
        double specular, 
        double reflection,
        double shininess
    ) :
        r(r), g(g), b(b),
        ambient(ambient),
        diffuse(diffuse),
        specular(specular),
        reflection(reflection),
        shininess(shininess)
    { 
        // assert(isEqual(ambient + diffuse + specular + reflection), 1);
    }
};

struct Color {
    double r, g, b;

    Color() : r(0), g(0), b(0) { };

    Color(double r, double g, double b) :
        r(max(0.0, min(1.0, r))),
        g(max(0.0, min(1.0, g))),
        b(max(0.0, min(1.0, b)))    
    { }

    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }

    Color operator*(const Color& other) const {
        return Color(r * other.r, g * other.g, b * other.b);
    }

    Color operator+(const Color& other) const {
        return Color(r + other.r, g + other.g, b + other.b);
    }

    friend ostream& operator<<(ostream& out, const Color& color) {
        out << fixed << setprecision(2) << color.r << ", " << color.g << ", " << color.b;
        return out;
    }
};

// Global variables
double globalScale = 1.0;
double rotationAngle = 0.0;
bool showAxes = true;
double axesLength = 100.0;

double nearDistance, farDistance, fovY, aspectRatio;
int imageWidth, imageHeight;
int levelOfRecursion;

bool texture = false;
bool multiple = false;

struct Object;
struct Light;
struct Spotlight;

vector<unique_ptr<Object>> objects;
vector<Light> lights;
vector<Spotlight> spotlights;

vector<vector<Point>> pointBuffer;

struct Light {
    Point pos;
    double falloff;
    Color color;

    Light(const Point& pos, double falloff) :
        pos(pos), falloff(falloff), color(Color(1, 1, 1))
    {
        // cerr << "Created Light:\n";
        // cerr << "\tPos: " << pos << "\n";
        // cerr << "\tFalloff: " << falloff << "\n";
        // cerr << endl;
    }

    ~Light() {
        // cerr << "\nDestroying Light.\n" << endl;
    }
};

struct Spotlight {
    Point pos;
    double falloff;
    Point lookingAt;
    double cutoff;
    Color color;

    Spotlight(const Point& pos, double falloff, const Point& lookingAt, double cutoff) :
        pos(pos), falloff(falloff), lookingAt(lookingAt), cutoff(cutoff), color(Color(1, 1, 1))
    {
        // cerr << "Created Spotlight:\n";
        // cerr << "\tPos: " << pos << "\n";
        // cerr << "\tFalloff: " << falloff << "\n";
        // cerr << "\tLooking At: " << lookingAt << "\n";
        // cerr << "\tCutoff: " << cutoff << "\n";
        // cerr << endl;
    }

    ~Spotlight() {
        // cerr << "\nDestroying Spotlight.\n" << endl;
    }
};

struct Object {

    virtual ~Object() { }

    virtual void printObject() const { };

    virtual void SetColorProp(ColorProp cp) = 0;

    virtual ColorProp getColorProp() const = 0;

    virtual Color getColor(const Point& point) const = 0;

    virtual void draw() const = 0;

    virtual Vec3 getNormal(const Point& point, const Ray& ray) const = 0;

    virtual double getIntersection(const Ray& ray) const = 0;

    virtual Color getShading(
        const Point& intersectionPoint, 
        const Ray& ray,
        int levelOfRecursion
    ) const 
    {
        if (levelOfRecursion == 0) {
            return Color(0, 0, 0);
        }

        const auto color = getColor(intersectionPoint);
        const auto colorProp = getColorProp();

        const auto normal = getNormal(intersectionPoint, ray);
        auto reflectedRay = ray.getReflection(intersectionPoint, normal);
        reflectedRay.src = reflectedRay.src + reflectedRay.dir * EPSILON;

        auto ambientLight = color * colorProp.ambient;
        Color diffuseLight(0, 0, 0);
        Color specularLight(0, 0, 0);

        for (const Light& light : lights) {
            Vec3 toLightRayDir = light.pos - intersectionPoint;
            toLightRayDir.normalize();
            Point toLightRaySrc = intersectionPoint + toLightRayDir * EPSILON;
            Ray toLightRay(toLightRaySrc, toLightRayDir);

            double distanceToLight = (light.pos - toLightRaySrc).magnitude();
            bool free = true;
            for (int i = 0; i < objects.size(); ++i) {
                double t = objects[i]->getIntersection(toLightRay);
                if (t > 0 && t < distanceToLight) {
                    free = false;
                    break;
                }
            }
            if (!free) {
                continue;
            }

            double scalingFactor = exp(-(distanceToLight * distanceToLight * light.falloff));
            double lambert = max(0.0, toLightRayDir.dot(normal) * scalingFactor);
            double phong = max(0.0, pow(max(0.0, reflectedRay.dir.dot(toLightRay.dir)), colorProp.shininess) * scalingFactor);

            diffuseLight = diffuseLight + color * (colorProp.diffuse * lambert);
            specularLight = specularLight + light.color * (colorProp.specular * phong);
        }

        for (const Spotlight& spotlight : spotlights) {
            Vec3 toLightRayDir = spotlight.pos - intersectionPoint;
            toLightRayDir.normalize();
            Point toLightRaySrc = intersectionPoint + toLightRayDir * EPSILON;
            Ray toLightRay(toLightRaySrc, toLightRayDir);

            Vec3 lightDirection = spotlight.lookingAt - spotlight.pos;
            lightDirection.normalize();
            Vec3 lightRayDir = intersectionPoint - spotlight.pos;
            lightRayDir.normalize();
            double angle = RadToDeg(acos(lightRayDir.dot(lightDirection)));
            if (angle > spotlight.cutoff) {
                continue;
            }

            double distanceToLight = (spotlight.pos - toLightRaySrc).magnitude();
            bool free = true;
            for (int i = 0; i < objects.size(); ++i) {
                double t = objects[i]->getIntersection(toLightRay);
                if (t > 0 && t < distanceToLight) {
                    free = false;
                    break;
                }
            }
            if (!free) {
                continue;
            }

            double scalingFactor = exp(-(distanceToLight * distanceToLight * spotlight.falloff));
            double lambert = max(0.0, toLightRayDir.dot(normal) * scalingFactor);
            double phong = max(0.0, pow(max(0.0, reflectedRay.dir.dot(toLightRay.dir)), colorProp.shininess) * scalingFactor);

            diffuseLight = diffuseLight + color * (colorProp.diffuse * lambert);
            specularLight = specularLight + spotlight.color * (colorProp.specular * phong);
        }

        double tMin = -1;
        int objectIndex = -1;
        for (int k = 0; k < objects.size(); ++k) {
            double t = objects[k]->getIntersection(reflectedRay);

            if (t > 0 && (tMin == -1 || t < tMin)) {
                tMin = t;
                objectIndex = k;
            }
        }

        Color reflectedLight(0, 0, 0);
        if (tMin != -1) {
            Point reflectedIntersectinPoint = reflectedRay.src + reflectedRay.dir * tMin;
            reflectedLight = objects[objectIndex]->getShading(
                reflectedIntersectinPoint, 
                reflectedRay, 
                levelOfRecursion - 1
            );
        }
        reflectedLight = reflectedLight * colorProp.reflection; 

        auto shading = ambientLight + diffuseLight + specularLight + reflectedLight;

        return shading;
    }
};

struct Sphere : public Object {
    Point center;
    double radius;
    ColorProp colorProp;

    Sphere(const Point& center, double radius) : center(center), radius(radius) {
        // cerr << "\nCreating Sphere:\n";
        // cerr << "\tCenter: " << center << "\n";
        // cerr << "\tRadius: " << radius << "\n";
        // cerr << endl;
    }

    ~Sphere() {
        // cerr << "\nDestroying Sphere.\n" << endl;
    }

    void printObject() const {
        cerr << "Sphere:\n";
        cerr << "\tCenter: " << center << "\n";
        cerr << "\tRadius: " << radius << "\n";
        cerr << endl;
    }

    void SetColorProp(ColorProp cp) {
        double t = cp.ambient + cp.diffuse + cp.specular + cp.reflection;
        cp.ambient /= t;
        cp.diffuse /= t;
        cp.specular /= t;
        cp.reflection /= t;
        this->colorProp = cp;
    }

    ColorProp getColorProp() const {
        return colorProp;
    }

    Color getColor(const Point& point) const {
        return Color(colorProp.r, colorProp.g, colorProp.b);
    }

    void draw() const {
        glPushMatrix(); {
            glColor3d(colorProp.r, colorProp.g, colorProp.b);
            glTranslated(center.x, center.y, center.z);
            GLUquadric* quadric = gluNewQuadric();
            gluSphere(quadric, radius, 100, 20);
            gluDeleteQuadric(quadric);
        } glPopMatrix();
    }

    Vec3 getNormal(const Point& point, const Ray& ray) const {
        auto normal = point - center;
        normal.normalize();
        if (normal.dot(ray.dir) > 0) {
            normal = normal * -1;
        }
        return normal;
    }

    double getIntersection(const Ray& ray) const {
        // We need to transform the system so that center is at origin
        // i.e. translate by -center

        Vec3 raySrc = ray.src - center; // Translate ray source
        Vec3 rayDir = ray.dir; // Remains the same after translation

        double a = 1;
        double b = 2 * rayDir.dot(raySrc);
        double c = raySrc.dot(raySrc) - radius * radius;

        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return -1;
        }

        double t1 = (-b - sqrt(discriminant)) / (2 * a);
        double t2 = (-b + sqrt(discriminant)) / (2 * a);

        double t = -1;
        if (t1 > 0 && t2 > 0) {
            t =  min(t1, t2);
        } 
        else if (t1 > 0) {
            t = t1;
        } 
        else if (t2 > 0) {
            t =  t2;
        }

        if (t == -1) {
            return -1;
        }

        return t;
    }
};

struct Triangle : public Object {
    Point a, b, c;
    ColorProp colorProp;

    Triangle(const Point& a, const Point& b, const Point& c ) :
        a(a), b(b), c(c) 
    {
        // cerr << "\nCreating Triangle:\n";
        // cerr << "\ta: " << a << "\n";
        // cerr << "\tb: " << b << "\n";
        // cerr << "\tc: " << c << "\n";
        // cerr << endl;
    }

    ~Triangle() {
        // cerr << "\nDestroying Triangle.\n" << endl;
    }

    void printObject() const {
        // cerr << "Triangle:\n";
        // cerr << "\ta: " << a << "\n";
        // cerr << "\tb: " << b << "\n";
        // cerr << "\tc: " << c << "\n";
        // cerr << endl;
    }

    void SetColorProp(ColorProp cp) {
        double t = cp.ambient + cp.diffuse + cp.specular + cp.reflection;
        cp.ambient /= t;
        cp.diffuse /= t;
        cp.specular /= t;
        cp.reflection /= t;
        this->colorProp = cp;
    }

    ColorProp getColorProp() const {
        return colorProp;
    }

    Color getColor(const Point& point) const {
        return Color(colorProp.r, colorProp.g, colorProp.b);
    }

    void draw() const {
        glColor3d(colorProp.r, colorProp.g, colorProp.b);
        glBegin(GL_TRIANGLES); {
            glVertex3d(a.x, a.y, a.z);
            glVertex3d(b.x, b.y, b.z);
            glVertex3d(c.x, c.y, c.z);
        } glEnd();
    }

    Vec3 getNormal(const Point& point, const Ray& ray) const {
        Vec3 u = b - a;
        Vec3 v = c - a;
        Vec3 normal = u.cross(v); 
        normal.normalize();
        if (normal.dot(ray.dir) > 0) {
            normal = normal * -1;
        }
        return normal;
    }

    double getIntersection(const Ray& ray) const {
        auto determinant = [](const double mat[3][3]) {
            return (
                mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) 
              - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) 
              + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0])
            );
        };

        double AMatrix[3][3] = {
            {a.x - b.x, a.x - c.x, ray.dir[0]},
            {a.y - b.y, a.y - c.y, ray.dir[1]},
            {a.z - b.z, a.z - c.z, ray.dir[2]}
        };
        double betaMatrix[3][3] = {
            {a.x - ray.src.x, a.x - c.x, ray.dir[0]},
            {a.y - ray.src.y, a.y - c.y, ray.dir[1]},
            {a.z - ray.src.z, a.z - c.z, ray.dir[2]}
        };
        double gammaMatrix[3][3] = {
            {a.x - b.x, a.x - ray.src.x, ray.dir[0]},
            {a.y - b.y, a.y - ray.src.y, ray.dir[1]},
            {a.z - b.z, a.z - ray.src.z, ray.dir[2]}
        };
        double tMatrix[3][3] = {
            {a.x - b.x, a.x - c.x, a.x - ray.src.x},
            {a.y - b.y, a.y - c.y, a.y - ray.src.y},
            {a.z - b.z, a.z - c.z, a.z - ray.src.z}
        };

        double A = determinant(AMatrix);
        double beta = determinant(betaMatrix) / A;
        double gamma = determinant(gammaMatrix) / A;
        double t = determinant(tMatrix) / A;

        if (!(t > 0 && min(beta, gamma) > 0 && beta + gamma < 1)) {
            return -1;
        }

        return t;
    }
};

struct Checkerboard : public Object {
    double width;
    ColorProp colorProp;
    vector<vector<Color>> textureColor[2];

    Checkerboard(double width) :
        width(width)
    {
        // cerr << "\nCreating Checkerboard:\n";
        // cerr << "\tWidth = " << width << "\n";
        // cerr << endl;

        cerr << "\nReading texture_w.bmp" << endl;
        bitmap_image wImage("texture_w.bmp");
        textureColor[0].resize(wImage.height(), vector<Color>(wImage.width()));
        for (int x = 0; x < wImage.width(); ++x) {
            for (int y = 0; y < wImage.height(); ++y) {
                auto pixel = wImage.get_pixel(x, y);
                textureColor[0][y][x] = Color(pixel.red / 255.0, pixel.green / 255.0, pixel.blue / 255.0);
            }
        }
        cerr << "Finished reading texture_w.bmp\n" << endl;
        
        cerr << "\nReading texture_b.bmp" << endl;
        bitmap_image bImage("texture_b.bmp");
        textureColor[1].resize(bImage.height(), vector<Color>(bImage.width()));
        for (int x = 0; x < bImage.width(); ++x) {
            for (int y = 0; y < bImage.height(); ++y) {
                auto pixel = bImage.get_pixel(x, y);
                textureColor[1][y][x] = Color(pixel.red / 255.0, pixel.green / 255.0, pixel.blue / 255.0);
            }
        }
        cerr << "Finished reading texture_b.bmp\n" << endl;
    }

    ~Checkerboard() {
        // cerr << "\nDestroying Checkerboard.\n" << endl;
    }
    
    void printObject() const {
        cerr << "Checkerboard:\n";
        cerr << "\tWidth = " << width << "\n";
        cerr << endl;
    }

    void SetColorProp(ColorProp cp) {
        double t = cp.ambient + cp.diffuse + cp.specular + cp.reflection;
        cp.ambient /= t;
        cp.diffuse /= t;
        cp.specular /= t;
        cp.reflection /= t;
        this->colorProp = cp;
    }

    ColorProp getColorProp() const {
        return colorProp;
    }

    Color getColor(const Point& point) const {
        int r = floor(point.x / width);
        int c = floor(point.y / width);
        int k = ((r + c) % 2 + 2) % 2;
        if (texture) {
            return getTexture(point.x, point.y, k);
        }   
        if (k == 0) {
            return Color(0, 0, 0);
        }
        else {
            return Color(1, 1, 1);
        }
    }

    void draw() const {
        double cameraX = camera.pos.x, cameraY = camera.pos.y;
        int r = floor(cameraX / width);
        int c = floor(cameraY / width);
        int n = 32;
        for (int i = r - n; i < r + n; ++i) {
            for (int j = c - n; j < c + n; ++j) {
                Point points[4] = {
                    Point(i * width, j * width, 0),
                    Point((i + 1) * width, j * width, 0),
                    Point((i + 1) * width, (j + 1) * width, 0),
                    Point(i * width, (j + 1) * width, 0)                  
                };
                if ((i + j) % 2 == 0) {
                    glColor3d(0, 0, 0);
                }
                else {
                    glColor3d(1, 1, 1);
                }
                glBegin(GL_QUADS); {
                    for (int k = 0; k < 4; ++k) {
                        glVertex3d(points[k].x, points[k].y, points[k].z);
                    }
                } glEnd();
            }
        }
    }

    Vec3 getNormal(const Point& point, const Ray& ray) const {
        Vec3 normal({0, 0, 1});
        if (normal.dot(ray.dir) > 0) {
            normal = normal * -1;
        }
        return normal;
    }

    double getIntersection(const Ray& ray) const {
        auto normal = getNormal(Point(0, 0, 0), ray);

        double a = normal.dot(ray.dir);

        if (isEqual(a, 0)) {
            return -1;
        }

        double b = normal.dot(ray.src - Point(0, 0, 0));
        double t = -(b / a);

        if (t <= 0) {
            return -1;
        }

        return t;
    }

    Color getTexture(double x, double y, int t) const {
        int imageWidth = textureColor[t].front().size();
        int imageHeight = textureColor[t].size();
        double dx = width / imageWidth, dy = width / imageHeight;

        x = x - width * floor(x / width);
        y = y - width * floor(y / width);
        // assert(x >= 0 && y >= 0);

        int i = floor(x / dx);
        int j = floor(y / dy);
        j = imageHeight - 1 - j;

        return textureColor[t].at(j).at(i);
    }
};

// Draw axes
void drawAxes() {
    glLineWidth(10.0);
    glBegin(GL_LINES); {
        // X axis
        glColor3ub(255, 0, 0); // Red
        glVertex3f(-axesLength, 0, 0);
        glVertex3f(axesLength, 0, 0);

        // Y axis
        glColor3ub(0, 255, 0); // Green
        glVertex3f(0, -axesLength, 0);
        glVertex3f(0, axesLength, 0);

        // Z axis
        glColor3ub(0, 0, 255); // Blue
        glVertex3f(0, 0, -axesLength);
        glVertex3f(0, 0, axesLength);
    } glEnd();
}

// Read scene description.
void readDescription() {
    const string filename = "description.txt";
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file " + filename);
    }

    cerr << "\nReading from " << filename << "..." << endl;

    file >> nearDistance >> farDistance >> fovY >> aspectRatio;

    file >> levelOfRecursion;

    int screenSize;
    file >> screenSize;
    imageWidth = imageHeight = screenSize;

    double checkerboardSize;
    double checkerboardAmbient, checkerboardDiffuse, checkerboardReflection;
    file >> checkerboardSize;
    file >> checkerboardAmbient >> checkerboardDiffuse >> checkerboardReflection;
    auto checkerboard = make_unique<Checkerboard>(checkerboardSize);
    checkerboard->SetColorProp(ColorProp(
        1, 1, 1,
        checkerboardAmbient,
        checkerboardDiffuse,
        0,
        checkerboardReflection,
        0
    ));
    objects.push_back(move(checkerboard));

    int numberOfObjects;
    file >> numberOfObjects;

    for (int i = 0; i < numberOfObjects; ++i) {
        string type;
        file >> type;

        if (type == "sphere") {
            double centerX, centerY, centerZ;
            file >> centerX >> centerY >> centerZ;
            double radius;
            file >> radius;
            double colorRed, colorGreen, colorBlue;
            file >> colorRed >> colorGreen >> colorBlue;
            double ambient, diffuse, specular, reflection;
            file >> ambient >> diffuse >> specular >> reflection;
            double shininess;
            file >> shininess;

            auto sphere = make_unique<Sphere>(Point(centerX, centerY, centerZ), radius);
            sphere->SetColorProp(ColorProp(
                colorRed, colorGreen, colorBlue,
                ambient, diffuse, specular, reflection, shininess
            ));
            objects.push_back(move(sphere));
        }
        else if (type == "pyramid") {
            double lowPointX, lowPointY, lowPointZ;
            file >> lowPointX >> lowPointY >> lowPointZ;
            double width, height;
            file >> width >> height;
            double colorRed, colorGreen, colorBlue;
            file >> colorRed >> colorGreen >> colorBlue;
            double ambient, diffuse, specular, reflection;
            file >> ambient >> diffuse >> specular >> reflection;
            double shininess;
            file >> shininess;

            vector<Point> points;
            Point basePoint(lowPointX, lowPointY, lowPointZ);
            points.push_back(basePoint);
            points.push_back(basePoint + Vec3({width, 0, 0}));
            points.push_back(basePoint + Vec3({width, width, 0}));
            points.push_back(basePoint + Vec3({0, width, 0}));
            points.push_back(basePoint + Vec3({width/2, width/2, height}));

            int triangles[6][3] = {
                {0, 1, 4},
                {1, 2, 4},
                {2, 3, 4},
                {0, 4, 3},
                {0, 2, 1},
                {0, 3, 2}
            };

            for (int i = 0; i < 6; ++i) {
                int a = triangles[i][0], b = triangles[i][1], c = triangles[i][2];
                auto triangle = make_unique<Triangle>(points[a], points[b], points[c]);
                triangle->SetColorProp(ColorProp(
                    colorRed, colorGreen, colorBlue,
                    ambient, diffuse, specular, reflection, shininess
                ));
                objects.push_back(move(triangle));
            }
        }
        else if (type == "cube") {
            double buttomLeftX, buttomLeftY, buttomLeftZ;
            file >> buttomLeftX >> buttomLeftY >> buttomLeftZ;
            double sideLength;
            file >> sideLength;
            double colorRed, colorGreen, colorBlue;
            file >> colorRed >> colorGreen >> colorBlue;
            double ambient, diffuse, specular, reflection;
            file >> ambient >> diffuse >> specular >> reflection;
            double shininess;
            file >> shininess;

            vector<Point> points;
            points.push_back(Point(buttomLeftX, buttomLeftY, buttomLeftZ));
            points.push_back(points.front() + Vec3({sideLength, 0, 0}));
            points.push_back(points.front() + Vec3({sideLength, sideLength, 0}));
            points.push_back(points.front() + Vec3({0, sideLength, 0}));
            points.push_back(points.front() + Vec3({0, 0, sideLength}));
            points.push_back(points.front() + Vec3({sideLength, 0, sideLength}));
            points.push_back(points.front() + Vec3({sideLength, sideLength, sideLength}));
            points.push_back(points.front() + Vec3({0, sideLength, sideLength}));
            
            int triangles[12][3] = {
                {0, 2, 1}, {0, 3, 2},
                {4, 5, 6}, {4, 6, 7},
                {0, 1, 5}, {0, 5, 4},
                {2, 3, 6}, {3, 7, 6},
                {0, 4, 3}, {3, 4, 7},
                {2, 6, 5}, {2, 5, 1}
            };

            for (int i = 0; i < 12; ++i) {
                int a = triangles[i][0], b = triangles[i][1], c = triangles[i][2];
                auto triangle = make_unique<Triangle>(points[a], points[b], points[c]);
                triangle->SetColorProp(ColorProp(
                    colorRed, colorGreen, colorBlue,
                    ambient, diffuse, specular, reflection, shininess
                ));
                objects.push_back(move(triangle));
            }
        }
    }

    int numberOfLights;
    file >> numberOfLights;

    for (int i = 0; i < numberOfLights; ++i) {
        double posX, posY, posZ;
        file >> posX >> posY >> posZ;
        double falloff;
        file >> falloff;

        lights.emplace_back(Point(posX, posY, posZ), falloff);
        // lights.back().color = Color(1, 0, 0);
    }

    int numberOfSpotlight;
    file >> numberOfSpotlight;

    for (int i = 0; i < numberOfSpotlight; ++i) {
        double posX, posY, posZ;
        file >> posX >> posY >> posZ;
        double falloff;
        file >> falloff;
        double lookingAtX, lookingAtY, lookingAtZ;
        file >> lookingAtX >> lookingAtY >> lookingAtZ;
        double cutoffAngle;
        file >> cutoffAngle;

        spotlights.emplace_back(
            Point(posX, posY, posZ), 
            falloff, 
            Point(lookingAtX, lookingAtY, lookingAtZ), 
            cutoffAngle
        );
        // spotlights.back().color = Color(1, 0, 1);
    }

    cerr << "Finished reading " << filename << "\n" << endl;

    file.close();
}

void loadPointBuffer() {
    double planeWidth = 2 * nearDistance * tan(degToRad((aspectRatio * fovY) / 2));
    double planeHeight = 2 * nearDistance * tan(degToRad(fovY / 2));
    double dx = planeWidth / imageWidth;
    double dy = planeHeight / imageHeight;

    Point midPoint = camera.pos + (camera.look * nearDistance);
    Point topLeftPoint = midPoint - (camera.right * (planeWidth / 2)) + (camera.up * (planeHeight / 2));
    Point refPoint = topLeftPoint + (camera.right * (dx / 2)) - (camera.up * (dy / 2));

    pointBuffer.clear();
    pointBuffer.resize(imageHeight, vector<Point>(imageWidth));
    for (int i = 0; i < imageHeight; ++i) {
        for (int j = 0; j < imageWidth; ++j) {
            pointBuffer[i][j] = refPoint + (camera.right * (dx * j)) - (camera.up * (dy * i));
        }
    }
}

int frameCnt = 0;

void renderImage() {
    if (multiple) {
        ++frameCnt;
    }

    cerr << "\nRendering frame #" << frameCnt << endl;

    auto image = bitmap_image(imageWidth, imageHeight);

    Point origin(camera.pos.x, camera.pos.y, camera.pos.z);

    int totalPixel = imageWidth * imageHeight;
    int completedPixel = 0;

    loadPointBuffer();
    for (int i = 0; i < imageHeight; ++i) {
        for (int j = 0; j < imageWidth; ++j) {
            ++completedPixel;

            auto raySrc = pointBuffer[i][j];
            auto rayDir = pointBuffer[i][j] - origin;
            rayDir.normalize();
            Ray ray = Ray(raySrc, rayDir);

            double tMin = -1;
            int objectIndex = -1;
            for (int k = 0; k < objects.size(); ++k) {
                double t = objects[k]->getIntersection(ray);

                if (t > 0 && (tMin == -1 || t < tMin)) {
                    tMin = t;
                    objectIndex = k;
                }
            }

            if (tMin == -1 || tMin > farDistance) {
                continue;
            }

            Point intersectionPoint = ray.src + ray.dir * tMin;
            auto shading = objects[objectIndex]->getShading(intersectionPoint, ray, levelOfRecursion);

            image.set_pixel(j, i, 255 * shading.r, 255 * shading.g, 255 * shading.b);
        }
        double completed = ((double) completedPixel / totalPixel) * 100;
        cerr << "\rRendering " << fixed << setprecision(2) << completed << " %";
    }

    image.save_image(string("frame_") + to_string(frameCnt) + ".bmp");

    cerr << "\nDone rendering frame #" << frameCnt << "\n" << endl;
}

// Handler for window-repaint event.
// Call back when the window first appears
// and whenever the window needs to be re-painted.
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear

    glMatrixMode(GL_MODELVIEW); // To operate on Model-View matrix
    glLoadIdentity(); // Reset the model-view matrix

    camera.setupCamera();

    // Perform object rotation
    glRotated(rotationAngle, 0, 1, 0);

    // Perform global scaling
    glScaled(globalScale, globalScale, globalScale);

    if (showAxes) {
        drawAxes();
    }

    for (const auto& object : objects) {
        object->draw();
    }

    glColor3d(1, 1, 1);
    for (const auto& light : lights) {
        glPushMatrix(); {
            glTranslated(light.pos.x, light.pos.y, light.pos.z);
            GLUquadric* quadric = gluNewQuadric();
            gluSphere(quadric, 3, 50, 10);
            gluDeleteQuadric(quadric);
        } glPopMatrix();
    }
    for (const auto& spotlight : spotlights) {
        glPushMatrix(); {
            glTranslated(spotlight.pos.x, spotlight.pos.y, spotlight.pos.z);
            GLUquadric* quadric = gluNewQuadric();
            gluCylinder(quadric, 3, 0, 5, 50, 10);
            gluDeleteQuadric(quadric);
        } glPopMatrix();
    }

    glutSwapBuffers(); // Render
}

void idle() {
    glutPostRedisplay();
}

// Handler for window re-size event.
// Called back when the window first appears 
// and whenever the window is re-sized with its new width and height.
void reshapeListener(int width, int height) {
    // Compute new aspect ratio
    if (height == 0) {
        height = 1; // To prevent divide by 0
    }
    double aspect = (double) width / (double) height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity(); // Reset the projection matrix

    // Enable perspective projection with fovy, aspect, zNear and zFar
    // gluPerspective(45.0, aspect, 0.1, 100.0);
    gluPerspective(fovY, aspectRatio, nearDistance, farDistance);
}

// Callback handler for normal-key event.
void keyboardListener(unsigned char key, [[maybe_unused]] int x, [[maybe_unused]] int y) {
    switch (key) {
    case ' ':
        texture = !texture;
        if (texture) {
            cerr << "\nTexture on\n" << endl;
        }
        else {
            cerr << "\nTexture off\n" << endl;
        }
        break;

    case '0':
        renderImage();
        break;

    case '1':
        // std::cerr << "Looking left." << std::endl;
        camera.rotateLeft(ROTATION_ANGLE_DELTA);
        break;

    case '2':
        // std::cerr << "Looking right." << std::endl;
        camera.rotateRight(ROTATION_ANGLE_DELTA);
        break;

    case '3':
        // std::cerr << "Looking up." << std::endl;
        camera.rotateUp(ROTATION_ANGLE_DELTA);
        break;

    case '4':
        // std::cerr << "Looking down." << std::endl;
        camera.rotateDown(ROTATION_ANGLE_DELTA);
        break;

    case '5':
        // std::cerr << "Tilting counter-clockwise." << std::endl;
        camera.tiltClockwise(ROTATION_ANGLE_DELTA);
        break;

    case '6':
        // std::cerr << "Tilting clockwise." << std::endl;
        camera.tiltCounterclockwise(ROTATION_ANGLE_DELTA);
        break;

    case 'x':
        showAxes = !showAxes;
        break;

    case 'f':
        multiple = true;
        break;

    default:
        return;
    }

    glutPostRedisplay(); // Post a paint request to activate display()
}

// Callback handler for special-key event.
void specialKeyListener(int key, [[maybe_unused]] int x, [[maybe_unused]] int y) {
    switch (key) {
    case GLUT_KEY_PAGE_UP:
        // std::cerr << "Moving up." << std::endl;
        camera.moveUp(CAMERA_DELTA);
        break;

    case GLUT_KEY_PAGE_DOWN:
        // std::cerr << "Moving down." << std::endl;
        camera.moveDown(CAMERA_DELTA);
        break;

    case GLUT_KEY_LEFT:
        // std::cerr << "Moving left." << std::endl;
        camera.moveLeft(CAMERA_DELTA);
        break;

    case GLUT_KEY_RIGHT:
        // std::cerr << "Moving right." << std::endl;
        camera.moveRight(CAMERA_DELTA);
        break;

    case GLUT_KEY_UP:
        // std::cerr << "Moving forward." << std::endl;
        camera.moveForward(CAMERA_DELTA);
        break;

    case GLUT_KEY_DOWN:
        // std::cerr << "Moving backward." << std::endl;
        camera.moveBackward(CAMERA_DELTA);
        break;

    default:
        return;
    }
    
    glutPostRedisplay(); // Post a paint request to activate display()
}

// Initialize OpenGL.
void initGL() {
    glClearColor(0, 0, 0, 1);   // Set background color
    glEnable(GL_DEPTH_TEST);    // Enable depth testing for z-culling
}

// GLUT entry-point.
int main(int argc, char** argv) {
    readDescription();                          // Read scene description.
    glutInit(&argc, argv);                      // Initialize GLUT.
    glutInitWindowSize(imageWidth, imageHeight);// Set the window's initial width & height.
    glutInitWindowPosition(50, 50);             // Position the window's initial top-left corner.
    glutInitDisplayMode(
        GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);   // Enable Depth, Double buffer, RGB color.
    glutCreateWindow("Ray Tracing");            // Create a window with the given title.
    glutDisplayFunc(display);                   // Register display callback handler for window re-paint.
    glutIdleFunc(idle);                      // Register idle callback handler.
    glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape.
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event.
    glutSpecialFunc(specialKeyListener);        // Register callback handler for special-key event.
    initGL();                                   // OpenGL initialization.
    glutMainLoop();                             // Enter the event-processing loop

    return 0;
}
