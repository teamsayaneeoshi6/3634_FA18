// adapted from https://www.purplealienplanet.com/node/23

/* A simple ray tracer */
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h> /* Needed for boolean datatype */
#include <math.h>
#include <cuda.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

/* Width and height of out image */
#define SCALE 1
#define WIDTH  (2048*SCALE)
#define HEIGHT (1440*SCALE)
#define DEPTH  (2048*SCALE)

#define BOXSIZE (2048*SCALE)

#define TRIANGLE 1
#define SPHERE   2
#define CONE     3
#define DISK     4
#define CYLINDER 5
#define RECTANGLE 6

#define p_eps 1e-6

#define p_Nsamples 1

// ratio of importance in sampling primary ray versus random rays
#define p_primaryWeight 2.f

#define p_intersectDelta 0.1f

#define p_shadowDelta 0.1f
#define p_projectDelta 1e-2

#define p_maxLevel 4
#define p_maxNrays (2<<p_maxLevel)
#define p_maxNcollisions 8
#define p_apertureRadius 20.f
#define NRANDOM 10000



typedef struct {
  // uses bitfields (https://www.geeksforgeeks.org/bit-fields-c/)
  unsigned int refractor: 1; // 1 = launches refraction ray
  unsigned int reflector: 1; // 1 = launches reflection ray
  unsigned int emitter:   1; // 1 = emits light but does not launch forward ray
}info_t;


/* The vector structure */
typedef struct{
  dfloat x,y,z;
}vector_t;

/* The sphere */
typedef struct{
  vector_t pos;
  dfloat  radius;
  vector_t velocity;
  vector_t newVelocity;
  vector_t force;
}sphere_t; 

/* The triangle */
typedef struct{
  vector_t vertices[3];
}triangle_t;

/* The rectangle */
typedef struct{
  vector_t center;
  vector_t axis[2];
  dfloat   length[2]; // extent of rectangle in axis directions 
}rectangle_t;


/* The cone */
typedef struct{
  vector_t vertex; // apex
  vector_t axis;   // axis vector from apex into cone
  dfloat radius;   // base radius
  dfloat height;   // height of apex above base
}cone_t;

/* The disk */
typedef struct{
  vector_t center;
  vector_t normal;
  dfloat   radius;
}disk_t;

/* The cylinder */
typedef struct{
  vector_t center; // center of base disk face
  vector_t axis;   // axis vector from apex into cone
  dfloat radius;   // base radius
  dfloat height;   // height of apex above base
}cylinder_t;


typedef struct{
  dfloat xmin, xmax;
  dfloat ymin, ymax;
  dfloat zmin, zmax;
  int imin, imax;
  int jmin, jmax;
  int kmin, kmax;
}bbox_t;


/* union of shapes */
typedef struct {
  int id;

  int type;
  
  union {
    sphere_t   sphere;
    triangle_t triangle;
    cone_t     cone;
    disk_t     disk;
    cylinder_t cylinder;
    rectangle_t rectangle;
  };

  int material;  

  bbox_t bbox;
  
}shape_t;

/* The ray */
typedef struct{
  vector_t start;
  vector_t dir;
  int level; // which level of the recursion launched this ray
  dfloat   coef;
}ray_t;


/* Colour */
typedef struct{
  dfloat red, green, blue;
}colour_t;

/* Material Definition */
typedef struct{
  colour_t diffuse;
  dfloat reflection;
  dfloat refraction; // transmission coefficient
  dfloat eta;
  info_t info;
}material_t;

/* Lightsource definition */
typedef struct{
  vector_t pos;
  colour_t intensity;
}light_t;

/* sensor configuration */
typedef struct{
  vector_t eyeX;  // model world coordinates of eye
  vector_t Idir; // screen horizontaal direction in model world (unit)
  vector_t Jdir; // screen vertical direction in model world (unit)
  vector_t normal;  // normal through the screen in model world (unit)
  dfloat   Ilength; // length of screen horizontal direction in model world
  dfloat   Jlength; // length of screen vertical direction in model world
  dfloat   offset;
  colour_t bg;      // background color

  dfloat   focalPlaneOffset; // normal distance between sensor and focus plane
  vector_t lensC;       // center of thin lens
  
}sensor_t;

typedef struct{
  int NI; // number of cells in x direction
  int NJ; // number of cells in y direction
  int NK; // number of cells in z direction
  dfloat xmin;
  dfloat xmax;
  dfloat ymin;
  dfloat ymax;
  dfloat zmin;
  dfloat zmax;
  dfloat dx;
  dfloat dy;
  dfloat dz;
  dfloat invdx;
  dfloat invdy;
  dfloat invdz;

  int      boxCount;
  int *boxContents;
  bbox_t  *bboxes;
  int     *boxStarts;
  int *c_boxContents;
  bbox_t  *c_bboxes;
  int     *c_boxStarts;

}grid_t;

void saveppm(char *filename, unsigned char *img, int width, int height);



__host__ __device__ vector_t vectorCreate(dfloat x, dfloat y, dfloat z);
__host__ __device__ vector_t vectorSub(const vector_t v1, const vector_t v2);
__host__ __device__ dfloat   vectorDot(const vector_t v1, const vector_t v2);
__host__ __device__ vector_t vectorCrossProduct(const vector_t v1, const vector_t v2);
__host__ __device__ vector_t vectorScale(const dfloat c, const vector_t v);
__host__ __device__ vector_t vectorAdd(const vector_t v1, const vector_t v2);
__host__ __device__ dfloat   vectorTripleProduct(const vector_t a, const vector_t b, const vector_t c);
__host__ __device__ vector_t vectorNormalize(const vector_t a);
__host__ __device__ dfloat   vectorNorm(const vector_t a);
__host__ __device__ vector_t vectorOrthogonalize(const vector_t a, const vector_t b);

typedef struct{

  int Nmaterials;
  material_t *materials;
  material_t *c_materials;

  int Nshapes;
  shape_t *shapes;
  shape_t *c_shapes;
  
  int Nlights;
  light_t *lights;
  light_t *c_lights;

  grid_t grid;

  
  
} scene_t;

scene_t *sceneSetup();

void render(const scene_t *scene,
	    const dfloat costheta,
	    const dfloat sintheta,
	    unsigned char *img);




void saveImage(char *filename, unsigned char *img, int W, int H);


__device__ bool intersectRayTriangle(const ray_t r, const triangle_t tri, dfloat *t);
__device__ bool intersectRayRectangle(const ray_t r, const rectangle_t rect, dfloat *t);
__device__ bool intersectRaySphere(const ray_t r, const sphere_t s, dfloat *t);
__device__ bool intersectRayCone(const ray_t r, const cone_t cone, dfloat *t);
__device__ bool intersectRayDisk(const ray_t r, const disk_t disk, dfloat *t);
__device__ bool intersectRayCylinder(const ray_t r, const cylinder_t cylinder, dfloat *t);
__device__ bool intersectPointGridCell(const grid_t grid,
			    const vector_t p,
			    const int cellI,
			    const int cellJ,
			    const int cellK);
__device__ unsigned int intersectRayBox(ray_t *r, const bbox_t bbox);

__device__ bool intersectRayShape(const ray_t r, const shape_t *s, dfloat *t);
__host__ __device__ bool solveQuadratic(const dfloat a, const dfloat b, const dfloat c, dfloat *x0, dfloat *x1);

__host__ __device__ int iclamp(dfloat x, dfloat xmin, dfloat xmax);
dfloat clamp(dfloat x, dfloat xmin, dfloat xmax);

__host__ __device__ vector_t shapeComputeNormal(const vector_t v, const shape_t s);

__device__ material_t shapeComputeMaterial(const int Nmaterials, const material_t *materials,
				const vector_t v, const shape_t s);

bbox_t createBoundingBoxTriangle(triangle_t triangle);
bbox_t createBoundingBoxSphere(sphere_t sphere);
bbox_t createBoundingBoxCylinder(cylinder_t cylinder);
bbox_t createBoundingBoxCone(cone_t cone);
bbox_t createBoundingBoxRectangle(rectangle_t rectangle);
bbox_t createBoundingBoxDisk(disk_t disk);
bbox_t createBoundingBoxShape(const grid_t grid, shape_t shape);

void gridCountShapesInCellsKernel(const grid_t grid, const int Nshapes, shape_t *shapes, int *counts);

__device__ colour_t gridTrace(const grid_t grid,
		   const int Nshapes,
		   const shape_t *shapes,
		   const int Nlights,
		   const light_t *lights,
		   const int Nmaterials,
		   const material_t *materials,
		   ray_t  r,
		   int    level,
		   dfloat coef,
		   colour_t bg);

dfloat projectPointRectangle(const vector_t p, const rectangle_t rect, vector_t *closest);
dfloat projectPointDisk(const vector_t p, const disk_t disk, vector_t *closest);
dfloat projectPointCylinder(const vector_t p, const cylinder_t cylinder, vector_t *closest);
dfloat projectPointCone(vector_t p, const cone_t cone, vector_t *closest);
dfloat projectPointTriangle(vector_t p, const triangle_t triangle, vector_t *closest);
dfloat projectPointShape(const vector_t p, const shape_t shape, vector_t *closest);
	    

void interpolateScene(const int NI, const int NJ, unsigned char *img);

__device__ vector_t sensorLocation(const int NI,
			const int NJ,
			const int I,
			const int J,
			const sensor_t sensor);

void gridPopulate(grid_t *grid, int Nshapes, shape_t *shapes);

dim3 blockSize(NI,NJ,NK);
dim3 blockNum(WIDTH,HEIGHT,1);

__global__ void renderKernel << blockSize, blockNum >> (const int NI,
		  const int NJ,
		  const scene_t scene,
		  const sensor_t sensor,
		  const dfloat costheta,
		  const dfloat sintheta,
		  unsigned char *c_img);

void readPlyModel(const char *fileName, int *Ntriangles, triangle_t **triangles);

void sphereCollisions(const grid_t *grid,
		      const dfloat dt,
		      const dfloat g,
		      const int Nshapes,
		      shape_t *shapes);

void sphereUpdates(grid_t *grid,
		   const dfloat dt,
		   const dfloat g,
		   const int Nshapes,
		   shape_t *shapes);
