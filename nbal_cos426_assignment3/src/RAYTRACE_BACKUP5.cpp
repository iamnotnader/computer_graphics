// Source file for raytracing code


// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"

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
   	//int count = 0;
   
      double t = -(ray->Start().Vector().Dot(plane->Normal()) + plane->D())
         		  /(ray->Vector().Dot(plane->Normal()));	
      //cout <<"T "<<t<<endl;
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
      	//count++;
         inter.hit = false;
         return inter;
      }
   	
   	//if (count == 3 || count == 0)
   	//{
      inter.hit = true;
      inter.t = t;
      inter.pos = P;
      inter.node = node;
      inter.normal = plane->Normal();
      return inter;
   	//}
   	//else
   	//{
   	//inter.hit = false;
   	//return inter;
   	//}
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
      //double tmax = INT_MAX;
      R3Vector normal;
   	
      inter.hit = false;
   
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
		
		if (ray->Start().Z() < zmax && ray->Start().Z() > zmin
		&&  ray->Start().Y() < ymax && ray->Start().Y() > ymin
		&&  ray->Start().X() < xmax && ray->Start().X() > xmin)
		{
			//cout << "inside"<<endl;
			inter.hit = true;
         inter.t = 0;
         inter.normal = R3Vector(0,0,0);
         inter.node = node;
         inter.pos = R3Point(ray->Start());
			return inter;
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
      //cout << "CONE" <<endl;
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
   	
   	
   	//cout <<xp0<<" "<< yp0<<" "<< zp0<<" "<< vx<<" "<< vy<<" "<< vz<<" "<< xc<<" "<< yc<<" "<< zc<<" "<< c<<" "<< determinant<<endl;
   	
      if (determinant < 0)
      {
         //cout << "1"<<endl;
         inter.hit = false;
         return inter;
      }
   	
      double t2 =   (vx*(xc - xp0) + pow(c,2)*vy*(-yc + yp0) + vz*(zc - zp0) + 
         		sqrt(determinant))/(pow(vx,2) - pow(c,2)*pow(vy,2) + pow(vz,2));
      double t1 =  (vx*(xc - xp0) + c*c*vy*(-yc + yp0) + vz*(zc - zp0) - 
         		sqrt(determinant))/(vx*vx - c*c*vy*vy + vz*vz);
   	
   	//cout << t1 << " " << t2 << endl;
   	
   	//		cout << t1 <<endl;
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
      //cout << "CYLINDER" <<endl;
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
   	
      //cout <<xp0 << " " << zp0 << " " << vx << " " << vz<<" "<<r<<" "<<xc<<" "<<zc<<endl;
   	
      if (determinant < 0)
      {
         //cout << "1"<<endl;
         inter.hit = false;
         return inter;
      }
   	
      double t1 = (1/(vx*vx + vz*vz)) * (vx * (xc - xp0) + vz*zc - vz*zp0 - sqrt(determinant));
      double t2 = (1/(vx*vx + vz*vz)) * (vx * (xc - xp0) + vz*zc - vz*zp0 + sqrt(determinant));
      //cout <<xp0 << " " << zp0 << " " << vx << " " << vz<<" "<<r<<" "<<xc<<" "<<zc<<endl;
   	
      //cout << t1 <<" "<<t2<< endl;
   	
      if (t1 > 0)
      {
         //cout << "2"<<endl;
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
         //R3Point temp = R3Point(xc, inter.pos.Y(), zc);
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
      
      
      /*else if (t1 < 0 && t2 > 0)
      {
         //cout << "3"<<endl;
         inter.pos = ray->Start() + ray->Vector()*t2;
      	
         if (abs(inter.pos.Y() - yc) > cyl->Height()/2)
         {
            R3Intersection baseLow;
            R3Intersection baseHigh;
            baseLow = intersectCircle(ray, cyl->Center(), -cyl->Height()/2, cyl->Radius());
            baseHigh = intersectCircle(ray, cyl->Center(),  cyl->Height()/2, cyl->Radius());
         	
            baseLow.normal = R3Vector(0, -1, 0);
            baseHigh.normal = R3Vector(0, 1, 0);
            baseLow.node = node;
            baseHigh.node = node;
         	
            if (baseLow.hit == true)
            {
               if (baseHigh.hit == true)
               {
                  if (baseLow.t < baseHigh.t)
                     return baseLow;
                  else
                     return baseHigh;
               }
               return baseLow;
            }
            
            if (baseHigh.hit == true)
            {
               if (baseLow.hit == true)
               {
                  if (baseHigh.t < baseLow.t)
                     return baseHigh;
                  else
                     return baseLow;
               }
               return baseHigh;
            }
         	
            inter.hit = false;
            return inter;
         }
      	
      	
         inter.t = t2;
         inter.node = node;
         R3Point temp = R3Point(xc, inter.pos.Y(), zc);
         inter.normal = -(inter.pos - temp);
         inter.normal.Normalize();
      }
      */
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
      //cout<<"IntersectSphere"<<endl;
   	//Declarations
      R3Intersection inter;
      R3Sphere* sph = node->shape->sphere;
   	
      //cout << ray->Vector().X() << " " << ray->Vector().Y()<<endl;
   	
   	//Vector from start of ray to center of circle
      R3Vector L = sph->Center() - ray->Start();
   	
   	//Initialize hit to true.
      inter.hit = true;
   	
   	//If distance from start to center is negative, no hit.
      double tca = L.Dot(ray->Vector());
      if (tca < 0)
      {
         //cout <<"2"<<endl;
         inter.hit = false;
         return inter;
      }
   	
   	//If the distance between ray and center of circle > r, no hit.
      double d2 = L.Dot(L) - pow(tca, 2);
      double r2 = pow(sph->Radius(), 2);
      if (d2 > r2)
      {
         //cout << "1"<<endl;
         inter.hit = false;
         return inter;
      }
   	
   	//Otherwise two values of t.
      double thc = pow(r2 - d2, .5);
      double t1 = tca - thc;
      //double t2 = tca + thc;
   	
   	//If closer t is positive, hit (behind circle)
      if (t1 > 0)
      {	
         inter.t = t1;
         inter.pos = ray->Start() + ray->Vector()*t1;
         inter.node = node;
         inter.normal = (inter.pos - sph->Center());
         inter.normal.Normalize();
      }
      //If closer t is neg and farther t is pos, hit (inside circle)
      /*else if (t1 < 0 && t2 > 0)
      {
         inter.t = t2;
         inter.pos = ray->Start() + ray->Vector()*t2;
         inter.node = node;
         inter.normal = -(inter.pos - sph->Center());
         inter.normal.Normalize();
      }*/
   	
   
   	////cout<<"IntersectSphere"<<endl;
   	//Return the intersection.
      return inter;
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
      	//cout << "sphere";
         inter = IntersectSphere(node, &transformedRay);
      }
      else if ((node->shape != NULL) && (node->shape->type == R3_BOX_SHAPE))
      {
      	//cout << "sphere";
         //cout <<"tr " << p1.X() << " "<<p1.Y() << " " <<p1.Z()<<endl;
         ///cout << "vec "<<vec.X() << " " <<vec.Y()<<" " <<vec.Z()<<endl; 
         inter = IntersectBox(node, &transformedRay);
      	//cout << "-"<<endl;
      }
      else if ((node->shape != NULL) && (node->shape->type == R3_MESH_SHAPE))
      {
      	//cout <<"tr " << p1.X() << " "<<p1.Y() << " " <<p1.Z()<<endl;
         //cout << "vec "<<vec.X() << " " <<vec.Y()<<" " <<vec.Z()<<endl;
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
   	//inter.hit = false;
   	//cout << node->children.size()<<endl;
   	
      if (node->children.size() == 0)
         return inter;
      else
      {
         R3Intersection temp;
      	/*--------------------------------------*/
         for (unsigned int i = 0; i < node->children.size(); i++)
         {
         	/*COMMENT OUT TO REMOVE BBOX CHECKING.*/
            R3Node bbox;
            R3Shape shape;
            shape.box = &(node->bbox);
         
            bbox.shape = &shape;
            
            temp = IntersectBox(&bbox, ray);
            
            if ((temp.hit == true) && (temp.t < inter.t))
            {
            /*-------------------------------------*/
               temp = IntersectScene(node->children[i], 
                     	&transformedRay);
            	
               if (temp.hit == true)
                  if (temp.t < inter.t)
                     inter = temp;
            }//REMOVE THIS AS WELL.
         }
      }
   	
   	//
      R3Point pos = inter.pos;
      p1 = node->transformation * inter.pos;
      p2 = inter.pos + inter.normal;
      inter.pos = p1;
      inter.t = (inter.pos - ray->Start()).Length();
      inter.normal = (node->transformation * p2) - (node->transformation * pos);
      inter.normal.Normalize();
   	//inter.hit = false;
   	//cout<<"done"<<endl;
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
      //cout << "GeneratePixelRay" <<endl;
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
   
      ray = R3Ray(from, ((toRight-from)+(toUp-from)));
   
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
      R2Pixel square = R2Pixel(0, 0, 0, 0);
      R3Material* mat;
      R3Light* light;
      R3Vector reflect;
   			
      if (inter.hit == true)
      {
         if (inter.node != NULL)
            mat = inter.node->material;
         square = R2Pixel(mat->ka * scn->ambient + mat->emission);
      
         for (int i = 0; i < scn->NLights(); i++)
         {
         //cout << i;
            light = scn->Light(i);
            if (light->type == R3_DIRECTIONAL_LIGHT)
            {
               R3Vector L = light->direction;
               L.Flip();
               if (inter.normal.Dot(L) > 0)
               {
               	
                  R3Ray shadowRay = R3Ray(inter.pos + inter.normal*1e-6, L);
                  R3Intersection temp = IntersectScene(scn->Root(), &shadowRay);   
                  
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
            else if (light->type == R3_POINT_LIGHT)
            {
               double ca = light->constant_attenuation;
               double la = light->linear_attenuation;
               double qa = light->quadratic_attenuation;
               R3Vector L = light->position - inter.pos;
               double d = L.Length();
               L.Normalize();
               
               if (inter.normal.Dot(L) > 0)
               {
               	
                  R3Ray shadowRay = R3Ray(inter.pos + inter.normal*1e-6, L);
                  R3Intersection temp = IntersectScene(scn->Root(), &shadowRay);            
     			
						//cout << temp.t <<" "<<d<<endl;          
					
                  if (temp.hit == false || (temp.t > d))
                  {
                  square += (mat->kd) * inter.normal.Dot(L) 
                        		  * (light->color) * 1.0 / (ca + la*d + qa*d*d);
                  
                  reflect = 2 * inter.normal.Dot(L)*inter.normal - L;
                  reflect.Normalize();
                  if ((L).Dot(reflect) > 0)
                     square += (mat->ks) * pow((-1*ray->Vector()).Dot(reflect), 
                           mat->shininess) * (light->color) * 1.0 / (ca + la*d + qa*d*d);
                  }
               }
            }
            else if (light->type == R3_SPOT_LIGHT)
            {
               R3Vector L = light->position - inter.pos;
               double d = L.Length();
               L.Normalize();
               
               if (acos(L.Dot((-1.0)*light->direction)) <= light->angle_cutoff)
               {
                  double ca = light->constant_attenuation;
                  double la = light->linear_attenuation;
                  double qa = light->quadratic_attenuation;
                  double sd = light->angle_attenuation;
               
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




   R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
   int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection)
   {
      //cout <<"RenderImage " <<endl; 
   
   // Allocate  image
      R2Image *image = new R2Image(width, height);
      if (!image) {
         fprintf(stderr, "Unable to allocate image\n");
         return NULL;
      }
   
      R3Intersection inter;
   // FILL IN YOUR CODE HERE
      for (int i = 0; i < width; i++)
         for (int j = 0; j < height; j++)
         {
            //cout <<scene->Camera().eye.X();
            R3Ray ray = GeneratePixelRay(&(scene->Camera()), i, j, width, height);
            inter = IntersectScene(scene->Root(), &ray);
            image->Pixels(i)[j] = GetColor(scene, &ray, inter);
         	//cout<<endl;
         }
   
   // Return image
      return image;
   }