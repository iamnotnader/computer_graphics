// Source file for raytracing code


// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include <iostream>
   using namespace std;



////////////////////////////////////////////////////////////////////////
// Create image from scene
//
// This is the main ray tracing function called from raypro
// 
// "width" and "height" indicate the size of the ray traced image
//   (keep these small during debugging to speed up your code development cycle)
//
// "max_depth" indicates the maximum number of secondary reflections/transmissions to trace for any ray
//   (i.e., stop tracing a ray if it has already been reflected max_depth times -- 
//   0 means direct illumination, 1 means one bounce, etc.)
//
// "num_primary_rays_per_pixel" indicates the number of random rays to generate within 
//   each pixel during antialiasing.  This argument can be ignored if antialiasing is not implemented.
//
// "num_distributed_rays_per_intersection" indicates the number of secondary rays to generate
//   for each surface intersection if distributed ray tracing is implemented.  
//   It can be ignored otherwise.
// 
////////////////////////////////////////////////////////////////////////

/*-----------------------------------------------------------------------------
 * Structure: R3Intersection
 * Scope: Local
 *	Points: 1
 
	Structure returned by a ray-intersection method. 
 *---------------------------------------------------------------------------*/
   struct R3Intersection
   {
      bool hit;
      R3Node* node;
      R3Point pos;
      R3Vector normal;
      double t;
   };

/*-----------------------------------------------------------------------------
 * Structure: inFace
 * Scope: Local
 *	Points: 1
 
 	Determine if the specified ray intersects the face with arbitrary number 
	of vertices.
 *---------------------------------------------------------------------------*/
   static R3Intersection inFace(vector<R3MeshVertex*> verts, 
   									 R3Plane* plane, R3Node* node, R3Ray* ray)
   {
      R3Intersection inter;
      R3Vector v2;
      inter.hit = false;
   
      double t = -(ray->Start().Vector().Dot(plane->Normal()) + plane->D())
         		  /(ray->Vector().Dot(plane->Normal()));	
  	  
		
      if (t < 0)
      {
         inter.hit = false;
         return inter;
      }
      R3Point P = ray->Start() + t * ray->Vector();
      for (unsigned int i = 0; i < verts.size() - 1; i++)
      {
         v2 = verts[i+1]->position - P;
         v2.Cross(verts[i]->position - P);
         v2.Normalize();
         if (ray->Vector().Dot(v2) < 0)
         {
            inter.hit = false;
            return inter;
         }
      }
   	
      v2 = verts[0]->position - P;
      v2.Cross(verts[verts.size() - 1]->position - P);
      v2.Normalize();
      if (ray->Vector().Dot(v2) < 0)
      {	
         inter.hit = false;
         return inter;
      }
   	
      inter.hit = true;
      inter.t = t;
      inter.pos = P;
      inter.node = node;
      inter.normal = plane->Normal();
      return inter;
   }


/*-----------------------------------------------------------------------------
 * Structure: IntersectMesh
 * Scope: Local
 *	Points: 1
 
 	Check to see if the ray intersects the mesh and if so return an 
	R3Intersection.
 *---------------------------------------------------------------------------*/
   static R3Intersection IntersectMesh(R3Node* node, R3Ray* ray)
   {
      R3Intersection inter;
      R3Intersection temp;  
   	
      inter.hit = false;
      inter.t = INT_MAX;
		
      for (unsigned int i = 0; i < node->shape->mesh->faces.size(); i++)
      {
         temp = inFace(node->shape->mesh->faces[i]->vertices, 
            			&(node->shape->mesh->faces[i]->plane),  
            		   node, ray);
         if (temp.hit == true)      
            if (temp.t < inter.t)
               inter = temp;
      }
      return inter;
   }
	
	
	
	
/*-----------------------------------------------------------------------------
 * Structure: inRect
 * Scope: Local
 *	Points: 1
 
	Determine if the point specified is inside the rectangle with bootom left
	at xmin ymin and top right at xmax ymax. 
	
	Note: x and y don't have to be on the x and y axis. The third value of, say,
	zmin and zmax, is assumed constant.
 *---------------------------------------------------------------------------*/
   static bool inRect(double xmin, double xmax, double ymin, double ymax, 
   						 double px, double py)
   {
      if ((px <= xmax) && (px >= xmin))
         if ((py <= ymax) && (py >= ymin))
            return true;
      return false;
   }

/*-----------------------------------------------------------------------------
 * Method: IntersectBox
 * Scope: Local
 *	Points: 1
 
	Return the intersection of the ray with the box contained in node.
 *---------------------------------------------------------------------------*/
   static R3Intersection IntersectBox(R3Node* node, R3Ray* ray)
   {
      R3Intersection inter;
      R3Box* box = node->shape->box;
      double xmin = box->XMin();
      double ymin = box->YMin();
      double zmin = box->ZMin();
      double xmax = box->XMax();
      double ymax = box->YMax();
      double zmax = box->ZMax();
      double t;
      R3Vector normal;
      inter.hit = false;
   
   
   
      if (ray->Start().Z() < zmax && ray->Start().Z() > zmin
      &&  ray->Start().Y() < ymax && ray->Start().Y() > ymin
      &&  ray->Start().X() < xmax && ray->Start().X() > xmin)
      {
         normal = R3Vector(0, 0, 1);
         if (ray->Vector().Dot(normal) < 0)
         {
            t = (zmin - ray->Start().Z())/(ray->Vector().Z());
            if (t > 0)
            {
               if (inRect(xmin, xmax, ymin, ymax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Y() + ray->Vector().Y()*t))
               {
                  inter.hit = true;
                  inter.t = t;
                  inter.normal = normal;
                  inter.node = node;
                  inter.pos = R3Point(ray->Start() + t * ray->Vector());
                  return inter;   
               }
            }
         }
      
      
         normal = R3Vector(0, 0, -1);
         if (ray->Vector().Dot(normal) < 0)
         {
            t = (zmax - ray->Start().Z())/(ray->Vector().Z());
            if (t > 0)
            {
               if (inRect(xmin, xmax, ymin, ymax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Y() + ray->Vector().Y()*t))
               {
                  inter.hit = true;
                  inter.t = t;
                  inter.normal = normal;
                  inter.node = node;
                  inter.pos = R3Point(ray->Start() + t * ray->Vector());
                  return inter;   
               }
            
            }
         }
      
         normal = R3Vector(-1, 0, 0);
         if (ray->Vector().Dot(normal) < 0)
         {
            t = (xmax - ray->Start().X())/(ray->Vector().X());
            if (t > 0)
            {
               if (inRect(zmin, zmax, ymin, ymax, ray->Start().Z() + ray->Vector().Z()*t, ray->Start().Y() + ray->Vector().Y()*t))
               {
                  inter.hit = true;
                  inter.t = t;
                  inter.normal = normal;
                  inter.node = node;
                  inter.pos = R3Point(ray->Start() + t * ray->Vector());
                  return inter;   
               }
            }
         }
      
         normal = R3Vector(1, 0, 0);
         if (ray->Vector().Dot(normal) < 0)
         {
            t = (xmin - ray->Start().X())/(ray->Vector().X());
            if (t > 0)
            {
               if (inRect(zmin, zmax, ymin, ymax, ray->Start().Z() + ray->Vector().Z()*t, ray->Start().Y() + ray->Vector().Y()*t))
               {
                  inter.hit = true;
                  inter.t = t;
                  inter.normal = normal;
                  inter.node = node;
                  inter.pos = R3Point(ray->Start() + t * ray->Vector());
                  return inter;   
               }
            }
         }
      
         normal = R3Vector(0, 1, 0);
         if (ray->Vector().Dot(normal) < 0)
         {
            t = (ymin - ray->Start().Y())/(ray->Vector().Y());
            if (t > 0)
            {
               if (inRect(xmin, xmax, zmin, zmax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Z() + ray->Vector().Z()*t))
               {
                  inter.hit = true;
                  inter.t = t;
                  inter.normal = normal;
                  inter.node = node;
                  inter.pos = R3Point(ray->Start() + t * ray->Vector());
                  return inter;   
               }
            }
            
         }
      
         normal = R3Vector(0, -1, 0);
         if (ray->Vector().Dot(normal) < 0)
         {
            t = (ymax - ray->Start().Y())/(ray->Vector().Y());
            if (t > 0)
            {
               if (inRect(xmin, xmax, zmin, zmax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Z() + ray->Vector().Z()*t))
               {
                  inter.hit = true;
                  inter.t = t;
                  inter.normal = normal;
                  inter.node = node;
                  inter.pos = R3Point(ray->Start() + t * ray->Vector());
                  return inter;
               }
            }     
         }
         
         inter.hit = false;
         return inter;
      }
   
   
   
   
   
   
   
      normal = R3Vector(0, 0, 1);
      if (ray->Vector().Dot(normal) < 0)
      {
         t = (zmax - ray->Start().Z())/(ray->Vector().Z());
         if (t > 0)
         {
            if (inRect(xmin, xmax, ymin, ymax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Y() + ray->Vector().Y()*t))
            {
               inter.hit = true;
               inter.t = t;
               inter.normal = normal;
               inter.node = node;
               inter.pos = R3Point(ray->Start() + t * ray->Vector());
               return inter;   
            }
         }
      }
   
   
      normal = R3Vector(0, 0, -1);
      if (ray->Vector().Dot(normal) < 0)
      {
         t = (zmin - ray->Start().Z())/(ray->Vector().Z());
         if (t > 0)
         {
            if (inRect(xmin, xmax, ymin, ymax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Y() + ray->Vector().Y()*t))
            {
               inter.hit = true;
               inter.t = t;
               inter.normal = normal;
               inter.node = node;
               inter.pos = R3Point(ray->Start() + t * ray->Vector());
               return inter;   
            }
            
         }
      }
      
      normal = R3Vector(-1, 0, 0);
      if (ray->Vector().Dot(normal) < 0)
      {
         t = (xmin - ray->Start().X())/(ray->Vector().X());
         if (t > 0)
         {
            if (inRect(zmin, zmax, ymin, ymax, ray->Start().Z() + ray->Vector().Z()*t, ray->Start().Y() + ray->Vector().Y()*t))
            {
               inter.hit = true;
               inter.t = t;
               inter.normal = normal;
               inter.node = node;
               inter.pos = R3Point(ray->Start() + t * ray->Vector());
               return inter;   
            }
         }
      }
      
      normal = R3Vector(1, 0, 0);
      if (ray->Vector().Dot(normal) < 0)
      {
         t = (xmax - ray->Start().X())/(ray->Vector().X());
         if (t > 0)
         {
            if (inRect(zmin, zmax, ymin, ymax, ray->Start().Z() + ray->Vector().Z()*t, ray->Start().Y() + ray->Vector().Y()*t))
            {
               inter.hit = true;
               inter.t = t;
               inter.normal = normal;
               inter.node = node;
               inter.pos = R3Point(ray->Start() + t * ray->Vector());
               return inter;   
            }
         }
      }
   	
      normal = R3Vector(0, 1, 0);
      if (ray->Vector().Dot(normal) < 0)
      {
         t = (ymax - ray->Start().Y())/(ray->Vector().Y());
         if (t > 0)
         {
            if (inRect(xmin, xmax, zmin, zmax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Z() + ray->Vector().Z()*t))
            {
               inter.hit = true;
               inter.t = t;
               inter.normal = normal;
               inter.node = node;
               inter.pos = R3Point(ray->Start() + t * ray->Vector());
               return inter;   
            }
         }
      }
      
      normal = R3Vector(0, -1, 0);
      if (ray->Vector().Dot(normal) < 0)
      {
         t = (ymin - ray->Start().Y())/(ray->Vector().Y());
         if (t > 0)
         {
            if (inRect(xmin, xmax, zmin, zmax, ray->Start().X() + ray->Vector().X()*t, ray->Start().Z() + ray->Vector().Z()*t))
            {
               inter.hit = true;
               inter.t = t;
               inter.normal = normal;
               inter.node = node;
               inter.pos = R3Point(ray->Start() + t * ray->Vector());
               return inter;
            }
         }     
      }
   	
         	
      return inter;
   }
   
/*-----------------------------------------------------------------------------
 * Method: IntersectCircle
 * Scope: Local
 *	Points: 1
 
	Given a ray, a center point, and a radius, tell if the ray intersects the
	circle.
 *---------------------------------------------------------------------------*/
   R3Intersection intersectCircle(R3Ray* ray, R3Point center, 
   				double offset, double rad)
   {
      R3Intersection inter;
      double xc = center.X();
      double yc = center.Y() + offset;
      double zc = center.Z();
   	
      double t = (yc - ray->Start().Y())/(ray->Vector().Y());
      double x = ray->Start().X() + ray->Vector().X()*t - xc;
      double z = ray->Start().Z() + ray->Vector().Z()*t - zc;
   	
      if ((x*x + z*z) > rad*rad)
      {
         inter.hit = false;
         return inter;
      }
   	
      inter.hit = true;
      inter.t = t;
      inter.pos = ray->Start() + ray->Vector()*t;
   	
      return inter;
   }

/*-----------------------------------------------------------------------------
 * Method: IntersectCone
 * Scope: Local
 *	Points: 1
 
	Return an R3Intersection structure if the ray intersects the cone.
 *---------------------------------------------------------------------------*/
   static R3Intersection IntersectCone(R3Node* node, R3Ray* ray)
   {
      R3Intersection inter;
      R3Cone* cone = node->shape->cone;
   
   	//Initialize hit to true.
      inter.hit = false;
      
      double xp0 = ray->Start().X(); 
      double yp0 = ray->Start().Y();
      double zp0 = ray->Start().Z();
      double vx = ray->Vector().X();
      double vy = ray->Vector().Y();
      double vz = ray->Vector().Z();
      double xc = cone->Center().X();
      double yc = cone->Center().Y() + .5*cone->Height();
      double zc = cone->Center().Z();
      double c = cone->Radius()/cone->Height();
      double determinant =  -pow(vz*(xc - xp0) + vx*(-zc + zp0),2) + 
         			pow(c,2)*((pow(vx,2) + pow(vz,2))*pow(yc - yp0,2) + 
         			pow(vy,2)*(pow(xc - xp0,2) + pow(zc - zp0,2)) - 
         			2*vy*(yc - yp0)*(vx*xc - vx*xp0 + vz*zc - vz*zp0));
   	
      if (determinant < 0)
      {
         inter.hit = false;
         return inter;
      }
   	
      double t2 =   (vx*(xc - xp0) + pow(c,2)*vy*(-yc + yp0) + vz*(zc - zp0) + 
         		sqrt(determinant))/(pow(vx,2) - pow(c,2)*pow(vy,2) + pow(vz,2));
      double t1 =  (vx*(xc - xp0) + c*c*vy*(-yc + yp0) + vz*(zc - zp0) - 
         		sqrt(determinant))/(vx*vx - c*c*vy*vy + vz*vz);
   	
      if (t1 > 0)
      {
         if ((ray->Start() + ray->Vector()*t1).Y() > yc && (ray->Start() 
         	+ ray->Vector()*t2).Y() < yc)
         {
            inter = intersectCircle(ray, cone->Center(), -cone->Height()/2, cone->Radius());
         	
            if (inter.hit == true)
            {
               inter.normal = R3Vector(0, -1, 0);
               inter.node = node;
               return inter;
            }
         
         }
      	
         inter.pos = ray->Start() + ray->Vector()*t1;
      
         if (inter.pos.Y() > yc)
         {
            inter.hit = false;
            return inter;
         }
      
         if (inter.pos.Y() < (yc - cone->Height()))
         {
            if ((ray->Start() + ray->Vector()*t1).Y() > (ray->Start() + ray->Vector()*t2).Y())
            {
               inter.hit = false;
               return inter;
            }
            inter = intersectCircle(ray, cone->Center(), -cone->Height()/2, cone->Radius());
         	
            if (inter.hit == true)
            {
               inter.normal = R3Vector(0, -1, 0);
               inter.node = node;
               return inter;
            }
            else
               return inter;
         }
      
         inter.hit = true;
         inter.t = t1;
         inter.node = node;
         R3Point temp = R3Point(xc, inter.pos.Y(), zc);
         R3Vector horizontal = (inter.pos - temp);
         horizontal.Normalize();
         horizontal = horizontal * cone->Height();
         R3Vector vertical = (cone->Center()+.5*cone->Height()*R3posy_vector - temp);
         vertical.Normalize();
         vertical = vertical * cone->Radius();
         inter.normal = vertical + horizontal;
         inter.normal.Normalize();
         return inter;
      	
      }
   	
      return inter;
   }


/*-----------------------------------------------------------------------------
 * Method: IntersectCylinder
 * Scope: Local
 *	Points: 1
 
	Return an R3Intersection structure if the ray intersects the sphere.
 *---------------------------------------------------------------------------*/
   static R3Intersection IntersectCylinder(R3Node* node, R3Ray* ray)
   {
      R3Intersection inter;
      R3Cylinder* cyl = node->shape->cylinder;
      
   	//Initialize hit to true.
      inter.hit = true;
      double xp0 = ray->Start().X(); 
      double zp0 = ray->Start().Z();
      double vx = ray->Vector().X();
      double vz = ray->Vector().Z();
      double r = cyl->Radius();
      double xc = cyl->Center().X();
      double yc = cyl->Center().Y();
      double zc = cyl->Center().Z();
      double determinant = r*r * (vx*vx + vz*vz) - pow(vz*xc - vz*xp0 - vx*zc + vx*zp0, 2);
   	
   	
      if (determinant < 0)
      {
         inter.hit = false;
         return inter;
      }
   	
      double t1 = (1/(vx*vx + vz*vz)) * (vx * (xc - xp0) + vz*zc - vz*zp0 - sqrt(determinant));
      double t2 = (1/(vx*vx + vz*vz)) * (vx * (xc - xp0) + vz*zc - vz*zp0 + sqrt(determinant));
   	
      if (t1 > 0)
      {
      inter.pos = ray->Start() + ray->Vector()*t1;
      	
         if (inter.pos.Y() > (yc+cyl->Height()/2))
         {
            inter = intersectCircle(ray, cyl->Center(),  cyl->Height()/2, cyl->Radius());
            if (inter.hit == true)
            {
               inter.normal = R3Vector(0, 1, 0);
               inter.node = node;
               return inter;
            }
            else
               return inter;
         }
         if (inter.pos.Y() < (yc-cyl->Height()/2))
         {
            inter = intersectCircle(ray, cyl->Center(), -cyl->Height()/2, cyl->Radius());
            if (inter.hit == true)
            {
               inter.normal = R3Vector(0, -1, 0);
               inter.node = node;
               return inter;
            }
            else
               return inter;
         }
      	
         inter.t = t1;
         inter.node = node;
         inter.normal = (inter.pos - cyl->Center());
         inter.normal.SetY(0);
         inter.normal.Normalize();
      }
      else if (t1 < 0 && t2 > 0)
      {
         if (ray->Start().Y() < (yc-cyl->Height()/2))
         {
            inter = intersectCircle(ray, cyl->Center(), -cyl->Height()/2, cyl->Radius());
            if (inter.hit == true)
            {
               inter.normal = R3Vector(0, -1, 0);
               inter.node = node;
               return inter;
            }
            else
               return inter;
         }
      	
         if (ray->Start().Y() > (yc+cyl->Height()/2))
         {
            inter = intersectCircle(ray, cyl->Center(), cyl->Height()/2, cyl->Radius());
            if (inter.hit == true)
            {
               inter.normal = R3Vector(0, 1, 0);
               inter.node = node;
               return inter;
            }
            else
               return inter;
         }
      }
      
      else
      {
         inter.hit = false;
         return inter;
      }
   	
      return inter;
   }
	
/*-----------------------------------------------------------------------------
 * Method: IntersectSphere
 * Scope: Local
 *	Points: 1
 
	Return an R3Intersection structure if the ray intersects the sphere.
 *---------------------------------------------------------------------------*/
   static R3Intersection IntersectSphere(R3Node* node, R3Ray* ray)
   {
   	//Declarations
      R3Intersection inter;
      R3Sphere* sph = node->shape->sphere;
   	
   	//Vector from start of ray to center of circle
      R3Vector L = sph->Center() - ray->Start();
   	
   	//Initialize hit to false.
      inter.hit = false;
   	
   	//If distance from start to center is negative, no hit.
      double tca = L.Dot(ray->Vector());
      if (tca < 0)
      {
         inter.hit = false;
         return inter;
      }
   	
   	//If the distance between ray and center of circle > r, no hit.
      double d2 = L.Dot(L) - pow(tca, 2);
      double r2 = pow(sph->Radius(), 2);
      if (d2 > r2)
      {
         inter.hit = false;
         return inter;
      }
   	
   	//Otherwise two values of t.
      double thc = pow(r2 - d2, .5);
      double t1 = tca - thc;
      double t2 = tca + thc;
   	
   	//If closer t is positive, hit (behind circle)
      if (t1 > 0)
      {	
         inter.hit = true;
         inter.t = t1;
         inter.pos = ray->Start() + ray->Vector()*t1;
         inter.node = node;
         inter.normal = (inter.pos - sph->Center());
         inter.normal.Normalize();
      }
      //If closer t is neg and farther t is pos, hit (inside circle)
      else if (t1 < 0 && t2 > 0)
      {
         inter.hit = true;
         inter.t = t2;
         inter.pos = ray->Start() + ray->Vector()*t2;
         inter.node = node;
         inter.normal = -(inter.pos - sph->Center());
         inter.normal.Normalize();
      }
   	
      return inter;
   }

/*
	Comparison function for sorting.
*/
   bool compare(pair<double, int> a, pair<double, int> b)
   {
      return a.first < b.first;
   }


/*-----------------------------------------------------------------------------
 * Method: IntersectScene
 * Scope: Local
 *	Points: 0
 
	Compute intersection of ray with a scene (recursive).
 *---------------------------------------------------------------------------*/
   static R3Intersection IntersectScene(R3Node* node, R3Ray* ray)
   {
   	//cout <<"IntersectScene"<<endl;
      struct R3Intersection inter;
      inter.hit = false;
      inter.t = INT_MAX;
   	//
      R3Point p1 = node->transformation.Inverse() * ray->Start();
      R3Point p2 = ray->Start() + ray->Vector();
      R3Vector vec = (node->transformation.Inverse() * p2) 
         	- (node->transformation.Inverse() * ray->Start());
   	
      R3Ray transformedRay = R3Ray(p1, vec);
   	
      if ((node->shape != NULL) && (node->shape->type == R3_SPHERE_SHAPE))
      {
         inter = IntersectSphere(node, &transformedRay);
      }
      else if ((node->shape != NULL) && (node->shape->type == R3_BOX_SHAPE))
      {
         inter = IntersectBox(node, &transformedRay);
      }
      else if ((node->shape != NULL) && (node->shape->type == R3_MESH_SHAPE))
      {
         inter = IntersectMesh(node, &transformedRay);
      }  
      else if ((node->shape != NULL) && (node->shape->type == R3_CYLINDER_SHAPE))
      {
         inter = IntersectCylinder(node, &transformedRay);
      }  
      else if ((node->shape != NULL) && (node->shape->type == R3_CONE_SHAPE))
      {
         inter = IntersectCone(node, &transformedRay);
      }  
   	
      if (node->children.size() == 0)
         return inter;
      else
      {
         R3Intersection temp;
         vector< pair<double, int> > tValues;
      	/*--------------------------------------*/
         for (unsigned int i = 0; i < node->children.size(); i++)
         {
         	/*COMMENT OUT TO REMOVE BBOX CHECKING.*/
            R3Node bbox;
            R3Shape shape;
            shape.box = &(node->bbox);
         
            bbox.shape = &shape;
            
            temp = IntersectBox(&bbox, ray);
            
            if (ray->Start().Z() < shape.box->ZMax() && ray->Start().Z() > shape.box->ZMin()
            &&  ray->Start().Y() < shape.box->YMax() && ray->Start().Y() > shape.box->YMin()
            &&  ray->Start().X() < shape.box->XMax() && ray->Start().X() > shape.box->XMin())
               temp.t = 0;	
					    	
            if (temp.hit == true)
            {
               pair<double,int> p = pair<double, int>(temp.t, i);
               tValues.push_back(p);
            }
         }
         sort(tValues.begin(), tValues.end(), compare);
      	
         for (unsigned int i = 0; i < tValues.size(); i++)
         {
            /*-------------------------------------*/
         	//change tValues[i].second to i
            temp = IntersectScene(node->children[tValues[i].second], 
                     	&transformedRay);
            if (temp.hit == true)
            {
               if (temp.t < inter.t)
                  inter = temp;
            		/*COMMENT THIS OUT AS WELL*/
               if ((i == tValues.size() - 1) 
               	||  (temp.t < tValues[i+1].first))
                  break;
            		/*-----------------------*/
            }
         } 
      }
   	
      R3Point pos = inter.pos;
      p1 = node->transformation * inter.pos;
      p2 = inter.pos + inter.normal;
      inter.pos = p1;
      inter.t = (inter.pos - ray->Start()).Length();
      inter.normal = (node->transformation * p2) - (node->transformation * pos);
      inter.normal.Normalize();
      return inter;
   }

/*-----------------------------------------------------------------------------
 * Method: GeneratePixelRay
 * Scope: Local
 *	Points: 1
 
	Generate a ray going from the eye in the camera through the point
	specified on the image.
 *---------------------------------------------------------------------------*/
   static R3Ray GeneratePixelRay(struct R3Camera* cam, int i, int j, int w, int h)
   {
    	//Some standard math that works.  
   	
      R3Ray ray;
      R3Point PLo;
      R3Point PHi;
   
   
      R3Point from = cam->eye;
      PLo = from + cam->towards - tan(2*cam->xfov)*cam->right;
      PHi = from + cam->towards + tan(2*cam->xfov)*cam->right;
      R3Point toRight = PLo + ((i + 0.5) / w) * (PHi - PLo);
   
      PLo = from + cam->towards - tan(2*cam->yfov)*cam->up;
      PHi = from + cam->towards + tan(2*cam->yfov)*cam->up;
      R3Point toUp = (PLo + ((j + 0.5) / h) * (PHi - PLo));
   
      R3Vector vec = (toRight-from)+(toUp-from);
      vec.Normalize();
   	
      ray = R3Ray(from, vec);
   
      return ray;
   }

/*-----------------------------------------------------------------------------
 * Method: GetColor
 * Scope: Local
 *	Points: 0
 
	Given a scene, a ray, and an intersection, determine the color of the 
	pixel.
 *---------------------------------------------------------------------------*/
   static R2Pixel GetColor(R3Scene* scn, R3Ray* ray, R3Intersection inter)
   {
   	//Declarations/Definitions
      R2Pixel square = R2Pixel(0, 0, 0, 0);
      R3Material* mat;
      R3Light* light;
      R3Vector reflect;
   			
      if (inter.hit == true)
      {
         if (inter.node != NULL)
            mat = inter.node->material;
         else
            return R2Pixel(0,0,0,0);
         square = R2Pixel(mat->ka * scn->ambient + mat->emission);
      
      
      	//Do a calculation for each light.
         for (int i = 0; i < scn->NLights(); i++)
         {
            light = scn->Light(i);
            
         	//Directional light calculation
            if (light->type == R3_DIRECTIONAL_LIGHT)
            {
               R3Vector L = light->direction;
               L.Flip();
               if (inter.normal.Dot(L) > 0)
               {
               	
               	//Cast a shadow ray
                  R3Ray shadowRay = R3Ray(inter.pos + inter.normal*1e-6, L);
                  R3Intersection temp = IntersectScene(scn->Root(), &shadowRay);   
                  
               	//If it didn't hit anything do the calculaion.
                  if (temp.hit == false)
                  {
                     square += (mat->kd) * inter.normal.Dot(L) 
                        		  * (light->color);
                  
                     reflect = 2 * inter.normal.Dot(L)*inter.normal - L;
                     reflect.Normalize();
                  
                     if ((-1.0*ray->Vector()).Dot(reflect) > 0)
                        square += (mat->ks) * pow((-1*ray->Vector()).Dot(reflect), 
                           mat->shininess) * (light->color);
                  }
               }
            }
            
            //Point light calculation.
            else if (light->type == R3_POINT_LIGHT)
            {
               double ca = light->constant_attenuation;
               double la = light->linear_attenuation;
               double qa = light->quadratic_attenuation;
               R3Vector L = light->position - inter.pos;
               double d = L.Length();
               L.Normalize();
               
            	//If it didn't hit anything do the calculaion.
               if (inter.normal.Dot(L) > 0)
               {
               	
                  R3Ray shadowRay = R3Ray(inter.pos + inter.normal*1e-6, L);
                  R3Intersection temp = IntersectScene(scn->Root(), &shadowRay);            
               
                  if (temp.hit == false || (temp.t > d))
                  {
                     square += (mat->kd) * inter.normal.Dot(L) 
                        		  * (light->color) * 1.0 / (ca + la*d + qa*d*d);
                  
                     reflect = 2 * inter.normal.Dot(L)*inter.normal - L;
                     reflect.Normalize();
                     if ((-1*ray->Vector()).Dot(reflect) > 0)
                        square += (mat->ks) * pow((-1*ray->Vector()).Dot(reflect), 
                           mat->shininess) * (light->color) * 1.0 / (ca + la*d + qa*d*d);
                  }
               }
            }
            
            //Spot light calculation.
            else if (light->type == R3_SPOT_LIGHT)
            {
               R3Vector L = light->position - inter.pos;
               double d = L.Length();
               L.Normalize();
               
            	//Check angle is within cutoff.
               if (acos(L.Dot((-1.0)*light->direction)) <= light->angle_cutoff)
               {
                  double ca = light->constant_attenuation;
                  double la = light->linear_attenuation;
                  double qa = light->quadratic_attenuation;
                  double sd = light->angle_attenuation;
               
               	//If the object is facing the light, do the calculation.
                  if ((inter.normal.Dot(L) > 0) && (L.Dot((-1.0)*light->direction) > 0))
                  {
                  	
                     R3Ray shadowRay = R3Ray(inter.pos + inter.normal*1e-6, L);
                     R3Intersection temp = IntersectScene(scn->Root(), &shadowRay);            
                  
                     if (temp.hit == false || (temp.t > d))
                     {
                        square += (mat->kd) * inter.normal.Dot(L) 
                           	  * (light->color) 
                           	  * pow(L.Dot((-1.0)*light->direction), sd) 
                           	  * 1.0 / (ca + la*d + qa*d*d);
                     
                        reflect = 2 * inter.normal.Dot(L)*inter.normal - L;
                        reflect.Normalize();
                        if (((ray->Vector().Dot(reflect)) > 0) && (L.Dot((-1)*light->direction) > 0))
                           square += (mat->ks) * pow((((-1)*ray->Vector()).Dot(reflect)), 
                              mat->shininess) * (light->color) 
                              * pow(L.Dot((-1)*light->direction), sd) 
                              * 1.0 / (ca + la*d + qa*d*d);
                     }
                  }
               }
            }
         }
         
      }
      else
         square = R2Pixel(0, 0, 0, 0);
      return square;
   }
   
/*-----------------------------------------------------------------------------
 * Method: GetSpecular
 * Scope: Local
 *	Points: 0
 
	Given a scene, a ray, and an intersection, a count, and a depth, find the 
	component of reflection for the object.
 *---------------------------------------------------------------------------*/
   static R2Pixel GetSpecular(R3Scene* scn, R3Ray* ray, R3Intersection inter, 
   						int count, int max_depth)
   {
   	//Base case
      if ((inter.hit == false)
      ||	 (count >= max_depth)
      ||  (inter.node == NULL)
      ||  (inter.node->material == NULL)
      ||  (inter.node->material->ks == R2Pixel(0,0,0,0))
      )
         return GetColor(scn, ray, inter);
      
   	//Generate the reflection ray
      R3Vector reflect = 2 * inter.normal.Dot(-1*ray->Vector())*inter.normal 
         						+ ray->Vector();
      R3Ray specularRay = R3Ray(inter.pos, reflect);
      
   	//Cast the ray into the scene
      R3Intersection specularRef = IntersectScene(scn->Root(), &specularRay);
      
   	//Cary the lighting calculation up.
      return GetColor(scn, ray, inter) + inter.node->material->ks 
                  * GetSpecular(scn, &specularRay, specularRef, count+1, max_depth);
   }
	
   
/*-----------------------------------------------------------------------------
 * Method: GetSpecular
 * Scope: Local
 *	Points: 0
 
	Given a scene, a ray, and an intersection, a count, and a depth, find the 
	component of reflection for the object.
 *---------------------------------------------------------------------------*/
   static R2Pixel GetPerfTrans(R3Scene* scn, R3Ray* ray, R3Intersection inter, 
   						int count, int max_depth, double ni)
   {
   	//If we don't hit something or cont = max_depth, just return the 
   	//specular color calculation.
      if (inter.hit == false
      ||	 count >= (max_depth)
      ||  inter.node == NULL
      ||  inter.node->material == NULL
      ||  inter.node->material->kt == R2Pixel(0,0,0,0))
         return GetSpecular(scn, ray, inter, 0, max_depth);
   		
      R3Ray throughRay;
      R3Intersection throughRef;
      
   	//Otherwise, if the material has an index of refraction specified,
   	//then cast a ray with an angle.
      if (inter.node->material->indexofrefraction != 1
      &&  inter.node->material->indexofrefraction != 0)
      {
      	//Some math to get the ray.
         double rindex = inter.node->material->indexofrefraction;
         double n = ni / rindex;
         R3Vector N = inter.normal;
         double cosI = -N.Dot(ray->Vector());
         double cosT2 = 1.0 - n * n * (1.0 - cosI * cosI);
         if (cosT2 > 0.0)
         {
         	//Creating the ray and intersecting it with the scene.
            R3Vector T = (n * ray->Vector()) + (n * cosI - sqrt(cosT2)) * N;
            throughRay = R3Ray(inter.pos + T*1e-6, T);
            throughRef = IntersectScene(scn->Root(), &throughRay);
         
         	//If the ray doesn't hit anything, then just return the 
         	//specular color calculation.
            if (throughRef.hit == false)
               return GetSpecular(scn, ray, throughRef, 0, max_depth);
         
         	//Cast the ray out of the object to complete the calculation.
         	
         	//Some math..
            rindex = throughRef.node->material->indexofrefraction;
            n = inter.node->material->indexofrefraction/rindex;
            N = throughRef.normal;
            cosI = -N.Dot(throughRay.Vector());
            cosT2 = 1.0 - n * n * (1.0 - cosI * cosI);
            if (cosT2 > 0.0)
            {
            	//Create the new ray and intersect it with the scene.
               T = (n * throughRay.Vector()) + (n * cosI - sqrt(cosT2)) * N;
               throughRay = R3Ray(inter.pos + T*1e-6, T);
               throughRef = IntersectScene(scn->Root(), &throughRay);
            
            	//Carry the lighting calculation of the new ray up the scene graph.
               return GetSpecular(scn, ray, inter, 0, max_depth) + inter.node->material->kt 
                  * GetPerfTrans(scn, &throughRay, throughRef, count+1, max_depth, 
                  inter.node->material->indexofrefraction);
            }
         }
      }
        
   	//If no index of refraction is specified, cast a ray through the object.
      if (count >= (max_depth-1))
         return GetSpecular(scn, ray, inter, 0, max_depth);
      throughRay = R3Ray(inter.pos+ray->Vector()*1e-6, ray->Vector());
      throughRef = IntersectScene(scn->Root(), &throughRay);
      throughRay = R3Ray(throughRef.pos + ray->Vector()*1e-6, ray->Vector());
      throughRef = IntersectScene(scn->Root(), &throughRay);
      	
      return GetSpecular(scn, ray, inter, 0, max_depth) + inter.node->material->kt 
                  * GetPerfTrans(scn, &throughRay, throughRef, count+2, max_depth, 1);
      
      return GetSpecular(scn, ray, inter, 0, max_depth);
   }
	

   R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
   int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection)
   {
   // Allocate  image
      R2Image *image = new R2Image(width, height);
      if (!image) {
         fprintf(stderr, "Unable to allocate image\n");
         return NULL;
      }
   
      R3Intersection inter;
   // FILL IN YOUR CODE HERE
      for (int i = 0; i < width; i++)
      {
         for (int j = 0; j < height; j++)
         {
            R3Ray ray = GeneratePixelRay(&(scene->Camera()), i, j, width, height);
            inter = IntersectScene(scene->Root(), &ray);
            if (inter.hit == true)
               image->Pixels(i)[j] = GetPerfTrans(scene, &ray, inter, 0, max_depth, 1);
               						 //GetSpecular(scene, &ray, inter, 0, max_depth);  
            else
               image->Pixels(i)[j] = R2Pixel(0, 0, 0, 0);
         }
         cout <<((int)((double)i/width*100))<<"\% ("<<i<<"/"<<width<<")"<<endl;
      }
   
   // Return image
      return image;
   }