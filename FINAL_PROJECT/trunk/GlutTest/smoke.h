// code adapted from www.videotutorialsrock.com

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include "GL/glut.h"
#endif

#include "imageloader.h"
#include "vec3f.h"
#include "R3/R3.h"


using namespace std;

const float PI = 3.1415926535f;

//Returns a random float from 0 to < 1
float randomFloat();

//Represents a single particle.
struct Particle {
	Vec3f pos;
	Vec3f velocity;
	Vec3f color;
	float timeAlive; //The amount of time that this particle has been alive.
	float lifespan;  //The total amount of time that this particle is to live.
	float size;
};

//Rotates the vector by the indicated number of degrees about the specified axis
Vec3f rotate(Vec3f v, Vec3f axis, float degrees);

//Returns the position of the particle, after rotating the camera
Vec3f adjParticlePos(Vec3f pos);

//Returns whether particle1 is in back of particle2
bool compareParticles(Particle* particle1, Particle* particle2);

const float GRAVITY = 3.0f;
const int NUM_PARTICLES = 100;
//The interval of time, in seconds, by which the particle engine periodically
//steps.
const float STEP_TIME = 0.01f;
//The length of the sides of the quadrilateral drawn for each particle.
const float PARTICLE_SIZE = 0.1f;


class ParticleEngine {
	private:
		GLuint textureId;
		Particle particles[NUM_PARTICLES];
		//The amount of time until the next call to step().
		float timeUntilNextStep;
		//The angle at which the fountain is shooting particles, in radians.
		float angle;
		
		Vec3f curColor() {
			Vec3f color;
			
			// gray
			color = Vec3f(0.745f, 0.745f, 0.745f);

			//Make sure each of the color's components range from 0 to 1
			for(int i = 0; i < 3; i++) {
				if (color[i] < 0) {
					color[i] = 0;
				}
				else if (color[i] > 1) {
					color[i] = 1;
				}
			}
			
			return color;
		}
		
		//Returns the average velocity of particles produced by the fountain.
		Vec3f curVelocity() {
			return Vec3f(2 * cos(PI/2), 2.0f, 2 * sin(PI));
		}
		
		//Alters p to be a particle newly produced by the fountain.
		void createParticle(Particle* p) {
			p->pos = Vec3f(0, 0, 0);
			p->velocity = curVelocity() + Vec3f(0.5f * randomFloat() - 0.25f,
												0.5f * randomFloat() - 0.25f,
												0.5f * randomFloat() - 0.25f);
			p->color = curColor();
			p->timeAlive = 0;
			p->lifespan = randomFloat() + 0.5;
			p->size = PARTICLE_SIZE;
		}
		
		//Advances the particle fountain by STEP_TIME seconds.
		void step() {			
			angle += 0.5f * STEP_TIME;
			while (angle > 2 * PI) {
				angle -= 2 * PI;
			}
			
			for(int i = 0; i < NUM_PARTICLES; i++) {
				Particle* p = particles + i;
				
				p->pos += p->velocity * STEP_TIME;
				p->velocity += Vec3f(-GRAVITY * STEP_TIME, 0.2*randomFloat()-0.1, 0.0f);
				p->timeAlive += STEP_TIME;
				if (p->timeAlive > p->lifespan) {
					createParticle(p);
				}
				p->size += 0.001 * pow(25,p->timeAlive);
			}
		}

	public:
		ParticleEngine(GLuint textureId1);
		
		//Advances the particle fountain by the specified amount of time.
		void advance(float dt);
		
		//Draws the particle fountain.
		void draw();
};


//Returns an array indicating pixel data for an RGBA image that is the same as
//image, but with an alpha channel indicated by the grayscale image alphaChannel
char* addAlphaChannel(Image* image, Image* alphaChannel);
//Makes the image into a texture, using the specified grayscale image as an
//alpha channel and returns the id of the texture
GLuint loadAlphaTexture(Image* image, Image* alphaChannel);
// draw the particles system
void drawParticles(ParticleEngine* input, double tx, double ty, double tz, double sx, double sy, double sz, R3Point *shipPos, double cullFront, double cullBack);
