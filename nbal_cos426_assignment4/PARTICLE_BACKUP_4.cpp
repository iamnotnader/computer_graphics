// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytracer.cpp"
#include <iostream>
   using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif

#define PI 3.14159



////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

   double RandomNumber(void)
   {
   #ifdef _WIN32
   // Seed random number generator
   static int first = 1;
   if (first) {
    srand(GetTickCount());
    first = 0;
   }
   
   // Return random number
   int r1 = rand();
   double r2 = ((double) rand()) / ((double) RAND_MAX);
   return (r1 + r2) / ((double) RAND_MAX);
   #else 
   // Seed random number generator
      static int first = 1;
      if (first) {
         struct timeval timevalue;
         gettimeofday(&timevalue, 0);
         srand48(timevalue.tv_usec);
         first = 0;
      }
   
   // Return random number
      return drand48();
   #endif
   }

   static int discrete(double* arr, int size) 
   {
        // precondition: sum of array entries equals 1
      double r = RandomNumber();
      double sum = 0.0;
      for (int i = 0; i < size; i++) 
      {
         sum += arr[i];
         if (sum >= r) 
            return i;
      }
      assert (false);
      return -1;
   }
	

   static R3Particle* generate_random_box_particle(R3ParticleSource* source)
   {
      R3Box* box = source->shape->box;
      int which_face = RandomNumber() * 6;
      R3Point pos;
      R3Vector dir;
   	
      if (which_face == 0)
      {
         pos = R3Point(box->XMin() + RandomNumber()*(box->XMax()-box->XMin()), 
            	 		  box->YMin() + RandomNumber()*(box->YMax()-box->YMin()),
            			  box->ZMin());
         dir = R3Vector(0, 0, -1);
      }
      else if (which_face == 1)
      {
         pos = R3Point(box->XMin() + RandomNumber()*(box->XMax()-box->XMin()), 
            	 		  box->YMin() + RandomNumber()*(box->YMax()-box->YMin()),
            			  box->ZMax());
         dir = R3Vector(0, 0, 1);
      }
      else if (which_face == 2)
      {
         pos = R3Point(box->XMin(), 
            	 		  box->YMin() + RandomNumber()*(box->YMax()-box->YMin()),
            			  box->ZMin() + RandomNumber()*(box->ZMax()-box->ZMin()));
         dir = R3Vector(-1, 0, 0);
      }
      else if (which_face == 3)
      {
         pos = R3Point(box->XMax(), 
            	 		  box->YMin() + RandomNumber()*(box->YMax()-box->YMin()),
            			  box->ZMin() + RandomNumber()*(box->ZMax()-box->ZMin()));
         dir = R3Vector(1, 0, 0);
      }
      else if (which_face == 4)
      {
         pos = R3Point(box->XMin() + RandomNumber()*(box->XMax()-box->XMin()), 
            	 		  box->YMin(),
            			  box->ZMin() + RandomNumber()*(box->ZMax()-box->ZMin()));
         dir = R3Vector(0, -1, 0);
      }
      else if (which_face == 5)
      {
         pos = R3Point(box->XMin() + RandomNumber()*(box->XMax()-box->XMin()), 
            	 		  box->YMax(),
            			  box->ZMin() + RandomNumber()*(box->ZMax()-box->ZMin()));
         dir = R3Vector(0, 1, 0);
      }
   					
      R3Particle* particle = new R3Particle();
      particle->position = pos;
      particle->velocity = source->velocity * dir;
      particle->mass = source->mass;
      particle->fixed = source->fixed;
      particle->drag = source->drag;
      particle->elasticity = source->elasticity;
      particle->lifetime = source->lifetime;
      particle->material = source->material; 
      
      return particle;
   }


////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////
   void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
   {
   // Generate new particles for every source
      R3ParticleSource* source;
      for (int i = 0; i < scene->NParticleSources(); i++)
      {
         source = scene->ParticleSource(i);
         int num_particles = source->elapsed_time * source->rate;
      	
         if ((num_particles) >= 1)
         {
            source->elapsed_time = 0;
         	
            for (int m = 0; m < num_particles; m++)
            {
               if (source->shape->type == R3_SEGMENT_SHAPE)
               {
                  R3Segment* seg = source->shape->segment;
                  R3Point pos = seg->Start() + RandomNumber() * (seg->End() - seg->Start());
                  R3Vector dir = seg->Vector();
                  dir.Cross(R3Vector(RandomNumber() - .5, RandomNumber() - .5, RandomNumber() - .5));
                  dir.Normalize();
               	
                  R3Particle* particle = new R3Particle();
                  particle->position = pos;
                  particle->velocity = source->velocity * dir;
                  particle->mass = source->mass;
                  particle->fixed = source->fixed;
                  particle->drag = source->drag;
                  particle->elasticity = source->elasticity;
                  particle->lifetime = source->lifetime;
                  particle->material = source->material; 
               
                  scene->particles.push_back(particle);
               }
            
               if (source->shape->type == R3_BOX_SHAPE)
               {
                  R3Particle* particle = generate_random_box_particle(source);
                  scene->particles.push_back(particle);
               }
               
               else if (source->shape->type == R3_MESH_SHAPE)
               {
                  R3Mesh* mesh = source->shape->mesh;
                  double area_sum = 0.0;
                  double face_prob[mesh->faces.size()];
               
                  for (unsigned int i = 0; i < mesh->faces.size(); i++)
                  {
                     face_prob[i] = mesh->faces[i]->Area();
                     area_sum += face_prob[i];
                  }
                  for (unsigned int i = 0; i < mesh->faces.size(); i++)
                     face_prob[i] /= area_sum;	
               	
                  int source_face = discrete(face_prob, mesh->faces.size());
                  int num_verts = mesh->faces[source_face]->vertices.size();
                
                  R3Point pos = R3Point(0,0,0);
                  R3Vector dir = R3Vector(0,0,0);
                  if (num_verts == 3)
                  {
                     double a = 0;
                     double b = 0;
                     double c;
                     while (a + b <= 1)
                     {
                        a = RandomNumber();
                        b = RandomNumber();  
                     }
                  
                     a = 1 - a;
                     b = 1 - b;
                     c = 1 - a - b;
                  
                     pos = a * mesh->faces[source_face]->vertices[0]->position
                        + b * mesh->faces[source_face]->vertices[1]->position
                        + c * mesh->faces[source_face]->vertices[2]->position;
                     mesh->faces[source_face]->vertices[0]->UpdateNormal();
                     mesh->faces[source_face]->vertices[1]->UpdateNormal();
                     mesh->faces[source_face]->vertices[2]->UpdateNormal();
                  
                     dir = a * mesh->faces[source_face]->vertices[0]->normal
                        + b * mesh->faces[source_face]->vertices[0]->normal
                        + c * mesh->faces[source_face]->vertices[0]->normal;
                  }
                  else
                  {
                     int num_verts = mesh->faces[source_face]->vertices.size();
                     double vert_sum = 0.0;
                     double vert_prob[num_verts];
                  
                     for (int i = 0; i < num_verts; i++)
                     {
                        vert_prob[i] = RandomNumber();
                        vert_sum += vert_prob[i];
                     }
                  
                     R3Point pos = R3Point(0,0,0);
                     R3Vector dir = R3Vector(0,0,0);
                     for (int i = 0; i < num_verts; i++)
                     {
                        vert_prob[i] /= vert_sum;
                        pos += vert_prob[i] * mesh->faces[source_face]->vertices[i]->position;
                        mesh->faces[source_face]->vertices[i]->UpdateNormal();
                        dir += vert_prob[i] * mesh->faces[source_face]->vertices[i]->normal;
                        //cout << vert_prob[i] << endl;
                     }
                  }
               
                  dir.Normalize();
               
                  R3Particle* particle = new R3Particle();
                  particle->position = pos;
                  particle->velocity = source->velocity * dir;
                  particle->mass = source->mass;
                  particle->fixed = source->fixed;
                  particle->drag = source->drag;
                  particle->elasticity = source->elasticity;
                  particle->lifetime = source->lifetime;
                  particle->material = source->material; 
               
                  scene->particles.push_back(particle);
               	
               }
               else if (source->shape->type == R3_SPHERE_SHAPE)
               {
                  R3Sphere* sphere = source->shape->sphere;
                  double radius = sphere->Radius();
                  R3Point center = sphere->Center();
               
                  double z = (RandomNumber() * 2 - 1);
                  double t = RandomNumber() * 2 * PI;
               
                  double r = sqrt(1.0 - z*z);
                  double x = r * cos(t) * radius;
                  double y = r * sin(t) * radius;
               
                  R3Point pos = R3Point(center.X() + x, center.Y() + y, center.Z() + z*radius);
                  R3Vector dir = (pos - center);
                  dir.Normalize();
               
                  R3Particle* particle = new R3Particle();
                  particle->position = pos;
                  particle->velocity = source->velocity * dir;
                  particle->mass = source->mass;
                  particle->fixed = source->fixed;
                  particle->drag = source->drag;
                  particle->elasticity = source->elasticity;
                  particle->lifetime = source->lifetime;
                  particle->material = source->material; 
               
                  scene->particles.push_back(particle);
               }
            }
         }
         else
            source->elapsed_time += delta_time;
      }
   }
	


   static R3Vector calculate_sink_force(R3Scene* scene, R3Particle* particle)
   {
      R3ParticleSink* sink;
      R3Vector force = R3Vector(0,0,0);  
   	
      for (int i = 0; i < scene->NParticleSinks(); i++)
      {
         sink = scene->ParticleSink(i);
      	
         if (sink->shape->type == R3_SPHERE_SHAPE)
         {
            R3Sphere* sphere = sink->shape->sphere;
            R3Point center = sphere->Center();
            
            R3Vector f = -(particle->position - center);
            double d = f.Length() - sphere->Radius();
            f.Normalize();
         	
            double mag = sink->intensity / (sink->constant_attenuation + 
               						sink->linear_attenuation*d + 
               						sink->quadratic_attenuation*d*d);
            force += f*mag;
         }
         else if (sink->shape->type == R3_MESH_SHAPE)
         {
            R3Mesh* mesh = sink->shape->mesh;
            R3Point center = R3Point(0,0,0);
            for (unsigned int j = 0; j < mesh->vertices.size(); j++)
            {
               center += mesh->vertices[i]->position / mesh->vertices.size();
            }
            R3Vector f = -(particle->position - center);
            double d = f.Length();
            f.Normalize();
         	
            double mag = sink->intensity / (sink->constant_attenuation + 
               						sink->linear_attenuation*d + 
               						sink->quadratic_attenuation*d*d);
            force += f*mag;
         }
         else if (sink->shape->type == R3_BOX_SHAPE)
         {
            R3Box* box = sink->shape->box;
            R3Point center = R3Point((box->XMax() - box->XMin())/2 + box->XMin(),
               							 (box->YMax() - box->YMin())/2 + box->YMin(),
               							 (box->ZMax() - box->ZMin())/2 + box->ZMin());
            R3Vector f = -(particle->position - center);
            double d = f.Length();
            f.Normalize();
         	
            double mag = sink->intensity / (sink->constant_attenuation + 
               						sink->linear_attenuation*d + 
               						sink->quadratic_attenuation*d*d);
            force += f*mag;
         
         }
      }
      return force;
   }
   
   static bool in_sink(R3Scene* scene, R3Particle* particle, double delta_time)
   {
      R3ParticleSink* sink;
   	
      for (int i = 0; i < scene->NParticleSinks(); i++)
      {
         sink = scene->ParticleSink(i);
      	
         if (sink->shape->type == R3_SPHERE_SHAPE)
         {
            R3Sphere* sphere = sink->shape->sphere;
            R3Point center = sphere->Center();
            double radius = sphere->Radius();  
         	
            if ((particle->position - center).Length() < radius)
               return true;
         }
         else if (sink->shape->type == R3_MESH_SHAPE)
         {
            R3Ray ray = R3Ray(particle->position, particle->velocity);
         	
            R3Node node;
            node.shape = sink->shape;
         	
            R3Intersection inter = IntersectMesh(&node, &ray);
         	
            if (inter.hit == true
            &&  inter.t <= ((particle->velocity).Length() * (delta_time) * 2))
               return true;
         }
         else if (sink->shape->type == R3_BOX_SHAPE)
         {
            R3Box* box = sink->shape->box;
            if (particle->position.Z() < box->ZMax() && particle->position.Z() > box->ZMin()
            &&  particle->position.Y() < box->YMax() && particle->position.Y() > box->YMin()
            &&  particle->position.X() < box->XMax() && particle->position.X() > box->XMin())
               return true;
         }
      }
      return false;
   }
	
   static R3Vector calculate_spring_force(R3Particle* p)
   {		
      R3Vector spring_force = R3Vector(0, 0, 0);
      R3Vector vp = p->velocity;
   	
      R3Vector D;
      double d;
      double ks;
      double kd;
      double s;
      R3Particle* q;
      R3Vector vq;
   	
      for (unsigned int i = 0; i < p->springs.size(); i++)
      {
         if (p->springs[i]->particles[0] == p)
            q = p->springs[i]->particles[1];
         else
            q = p->springs[i]->particles[0];
      	
         vq = q->old_velocity;
      	
         s = p->springs[i]->rest_length;
         ks = p->springs[i]->ks;
         kd = p->springs[i]->kd;
      	
         D = q->old_position - p->position;
         d = D.Length();
      	
         spring_force += (ks * (d - s) + kd * (vq - vp).Dot(D)) * D; 
      }
      
   	//cout << "spring_force: " << spring_force.X() << " " << spring_force.Y() << " " << spring_force.Z() << endl;
   	
      return spring_force;		
   }
   
   static R3Vector calculate_attractive_force(R3Scene* scene, R3Particle* particle)
   {
      R3Vector attractive_force = R3Vector(0,0,0);
   	
      double G = 6.67300e-11;
   	
      for (int i = 0; i < scene->NParticles(); i++)
      {
         if (scene->Particle(i) == particle)
            continue;
      	
         R3Vector r = particle->position - scene->Particle(i)->old_position;
         double m1 = particle->mass;
         double m2 = scene->Particle(i)->mass;
      	
         double mag = -G * m1 * m2 / (r.Length() * r.Length());
         r.Normalize();
      	
         attractive_force += mag * r;
      }
   	
      return attractive_force;
   }
	
	static R3Vector calculate_boid_force(R3Scene* scene, R3Particle* particle)
	{
		double range = 1;
		R3Vector net_force = R3Vector(0,0,0);
		vector<R3Particle*> neighbors;
		for (int i = 0; i < scene->NParticles(); i++)
		{
			if (scene->Particle(i) != particle
			&&  (scene->Particle(i)->old_position - particle->position).Length() <= range)
				neighbors.push_back(scene->Particle(i));
		}
		
		R3Point average_position = R3Point(0,0,0);
		R3Vector adjustment;
		
		for (unsigned int i = 0; i < neighbors.size(); i++)
			average_position += neighbors[i]->old_position/neighbors.size();
		adjustment = average_position - particle->position; 
		
		R3Vector toward_others = R3Vector(0,0,0);
		for (unsigned int i = 0; i < neighbors.size(); i++)
			toward_others += (neighbors[i]->old_position - particle->position)/pow((neighbors[i]->old_position - particle->position).Length(), 2); 
		
		R3Ray ray = R3Ray(particle->position, particle->velocity);
		R3Intersection inter = IntersectScene(scene->Root(), &ray);
		R3Vector avoid_collision = R3Vector(0,0,0);
		
		if (inter.hit == true)
			avoid_collision = 100*(inter.normal/(pow(inter.t, 2)+.1));
	
		return RandomNumber()*2*adjustment - toward_others + avoid_collision;
	}
	
   static R3Vector calculate_net_force(R3Scene* scene, R3Particle* particle)
   {
      R3Vector gravity = scene->gravity;
   
      R3Vector grav_force = particle->mass * gravity;
      R3Vector drag_force = -particle->drag * particle->velocity;
      R3Vector sink_force = calculate_sink_force(scene, particle);
      R3Vector spring_force = calculate_spring_force(particle);
      R3Vector attractive_force = calculate_attractive_force(scene, particle);//R3Vector(0,0,0);
   	R3Vector boid_force = calculate_boid_force(scene, particle);
   
      return grav_force + drag_force + sink_force + spring_force 
      	  + attractive_force + boid_force;
   }

   static void run_collisions(R3Scene* scene, R3Particle* particle, 
   							double* delta_time)
   {
      R3Ray ray = R3Ray(particle->position, particle->velocity);
      R3Intersection inter = IntersectScene(scene->Root(), &ray);
      
      while (inter.hit == true
      	&&  inter.t <= ((particle->velocity).Length() * (*delta_time) * 1.1))
      {
         R3Vector net_force = calculate_net_force(scene, particle);
      
         double time_until_collision = inter.t / (particle->velocity).Length();
         R3Vector parallel_component = (inter.normal).Dot(-particle->velocity) 
            								 * (inter.normal);
         R3Vector perp_component = -(-(particle->velocity) - parallel_component);
      	
      	
      	
         particle->velocity += (net_force / particle->mass) * time_until_collision * .99;
      	
         particle->velocity = parallel_component * particle->elasticity 
            			+ perp_component;
      	
         particle->position += (time_until_collision*.99) * particle->velocity; 	
      					
         *delta_time -= time_until_collision*.99;
         R3Ray r = R3Ray(particle->position, particle->velocity);
         inter = IntersectScene(scene->Root(), &r);
      }
   }


////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////
   void UpdateParticles(R3Scene *scene, double current_time, 
   						double delta_time, int integration_type)
   {
   // Update position for every particle
      R3Particle* particle;
      double error_threshold = 1e-15;  
   	
      R3Vector net_force;
   	
   	/**/
   	for (int h = 0; h < scene->NParticles(); h++)
         {
            scene->Particle(h)->old_velocity = scene->Particle(h)->velocity;
            scene->Particle(h)->old_position = scene->Particle(h)->position;
         }
   	
      particle = scene->Particle(RandomNumber() * scene->NParticles());
      net_force = calculate_net_force(scene, particle);
      double t_step = delta_time;
      R3Vector vel = particle->velocity;
      R3Point pos1 = R3Point(0, 0, 0);
      R3Point pos2 = R3Point(1, 1, 1);
      int div_num = 1;
      	
      while ((pos2 - pos1).Length() > error_threshold
         &&     div_num < 0)
      {
         pos1 = particle->position + t_step * vel;
         pos2 = particle->position + t_step/2 * vel + 
                   t_step/2 * (vel + t_step/2 * (net_force / particle->mass));
         div_num *= 2;
         t_step /= 2;		
      }
   	
      for (int k = 0; k < div_num; k++)
      {	
         for (int h = 0; h < scene->NParticles(); h++)
         {
            scene->Particle(h)->old_velocity = scene->Particle(h)->velocity;
            scene->Particle(h)->old_position = scene->Particle(h)->position;
         }
      	
      
         int i = 0;
         while (i < scene->NParticles())
         {
            particle = scene->Particle(i);
         
            if (particle->fixed)
            {
               particle->elapsed_time += delta_time;
               i++;
               continue;
            }
         	
         	//cout << time_passed << endl;
            if (in_sink(scene, particle, t_step)
            ||  particle->elapsed_time >= particle->lifetime
            &&  particle->lifetime != 0)
            {
               scene->particles.erase(scene->particles.begin() + i);
               i++;
               continue;
            }
         
            double time_remaining = t_step; 
         	
            run_collisions(scene, particle, &time_remaining);
         
            net_force = calculate_net_force(scene, particle);
         
            particle->position += time_remaining * particle->velocity;
            particle->velocity += time_remaining * (net_force / particle->mass);
         
            particle->elapsed_time += t_step;
            i++;
         }
      }
   }



////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

   void RenderParticles(R3Scene *scene, double current_time, double delta_time)
   {
   // Draw every particle
   
   // REPLACE CODE HERE
      glDisable(GL_LIGHTING);
      glPointSize(5);
      glBegin(GL_POINTS);
      for (int i = 0; i < scene->NParticles(); i++) {
         R3Particle *particle = scene->Particle(i);
         glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
         const R3Point& position = particle->position;
         glVertex3d(position[0], position[1], position[2]);
      }   
      glEnd();
   }



