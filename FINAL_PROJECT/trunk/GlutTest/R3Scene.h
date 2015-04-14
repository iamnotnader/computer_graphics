// Include file for the R3 scene stuff

#define R3Rgb R2Pixel

#include "R3/R3.h"

struct SFEnemy;
struct SFProjectile;

// Constant definitions

typedef enum {
    R3_BOX_SHAPE,
    R3_SPHERE_SHAPE,
    R3_CYLINDER_SHAPE,
    R3_CONE_SHAPE,
    R3_MESH_SHAPE,
    R3_SEGMENT_SHAPE,
    R3_NUM_SHAPE_TYPES
} R3ShapeType;

typedef enum {
    R3_DIRECTIONAL_LIGHT,
    R3_POINT_LIGHT,
    R3_SPOT_LIGHT,
    R3_AREA_LIGHT,
    R3_NUM_LIGHT_TYPES
} R3LightType;



// Scene element definitions

struct R3Shape {
    R3ShapeType type;
    R3Box *box;
    R3Sphere *sphere;
    R3Cylinder *cylinder;
    R3Cone *cone;
    R3Mesh *mesh;
    R3Segment *segment;
};  

struct R3Material {
    R3Rgb ka;
    R3Rgb kd;
    R3Rgb ks;
    R3Rgb kt;
    R3Rgb emission;
    double shininess;
    double indexofrefraction;
    R2Image *texture;
    int texture_index;
    int id;
};

struct R3Light {
    R3LightType type;
    R3Point position;
    R3Vector direction;
    double radius;
    R3Rgb color;
    double constant_attenuation;
    double linear_attenuation;
    double quadratic_attenuation;
    double angle_attenuation;
    double angle_cutoff;
};

struct R3Camera {
    R3Point eye;
    R3Vector towards;
    R3Vector right;
    R3Vector up;
    double xfov, yfov;
    double neardist, fardist;
};

struct R3Node {
    struct R3Node *parent;
    vector<struct R3Node *> children;
    R3Shape *shape;
    R3Matrix transformation;
    R3Matrix cumulativeTransformation; 
    // transformation including those of this node's parents; needed just for nodes with enemies
    R3Material *material;
    R3Box bbox;
    SFEnemy *enemy;
};


// Scene graph definition

struct R3Scene {
public:
    // Constructor functions
    R3Scene(void);
    
    // Access functions
    R3Node *Root(void) const;
    int NLights(void) const;
    R3Light *Light(int k) const;
    int NEnemies(void) const;
    SFEnemy *Enemy(int k) const;
    int NProjectiles(void) const;
    SFProjectile *Projectile(int k) const;
    R3Camera& Camera(void);
    R3Box& BBox(void);
    
    // I/O functions
    int Read(const char *filename, R3Node *root = NULL);
    
public:
    R3Node *root;
    R3Node *arwingNode;
    vector<R3Light *> lights;
    vector<SFEnemy *> enemies;
    vector<SFProjectile *> projectiles;
    R3Camera camera;
    R3Box bbox;
    R3Rgb background;
    R3Rgb ambient;
};



// Inline functions 

inline R3Node *R3Scene::
Root(void) const
{
    // Return root node
    return root;
}



inline int R3Scene::
NLights(void) const
{
    // Return number of lights
    return lights.size();
}



inline R3Light *R3Scene::
Light(int k) const
{
    // Return kth light
    return lights[k];
}


inline int R3Scene::
NEnemies(void) const
{
    // Return number of enemies
    return enemies.size();
}



inline SFEnemy *R3Scene::
Enemy(int k) const
{
    // Return kth enemy
    return enemies[k];
}

inline int R3Scene::
NProjectiles(void) const
{
    // Return number of projectiles
    return projectiles.size();
}



inline SFProjectile *R3Scene::
Projectile(int k) const
{
    // Return kth projectile
    return projectiles[k];
}


inline R3Camera& R3Scene::
Camera(void) 
{
    // Return camera
    return camera;
}



inline R3Box& R3Scene::
BBox(void) 
{
    // Return bounding box 
    return bbox;
}



//SFEnemy Declaration; I would much rather get it into a separate file
//but I can't get my head around these crazy linkage errors

struct SFEnemy {
    
    
public:
    //Constructors
    SFEnemy(void);
    SFEnemy(int fix, R3Mesh *mesh, R3Vector& initialVelocity, int health, float particle_velocity, int firing_rate);
    
    //  void Update(void);  //handled in gluttest.cpp
    
public:
    int fixed;              //does this enemy move?
    
    R3Point position;
    R3Vector movementPath;  //direction of enemy's movement, if applicable
    R3Vector towards;       //enemy's facing direction
    
    int firingRate;
    double movementSpeed;
    
    double fov;             //(radians) breadth of enemy's vision
    double firingSpread;    //(radians) breadth of enemy's firing range
    
    double projectileLength;
    double projectileSpeed;
    R3Point projectileSource;        //source of projectiles
    
    R3Node *node;
    R3Mesh *mesh;
    
    int health;
};

struct SFProjectile {
    //constructor
public:
    SFProjectile(void);
    SFProjectile(double spd, R3Node *pNode);
    
public:
    R3Segment segment;
    double power;
    double speed;
    R3Node *parentNode;
	 R3Node* node;
};
