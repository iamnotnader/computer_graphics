#include "smoke.h"


//Returns a random float from 0 to < 1
float randomFloat() {
	return (float)rand() / ((float)RAND_MAX + 1);
}

//Rotates the vector by the indicated number of degrees about the specified axis
Vec3f rotate(Vec3f v, Vec3f axis, float degrees) {
	axis = axis.normalize();
	float radians = degrees * PI / 180;
	float s = sin(radians);
	float c = cos(radians);
	return v * c + axis * axis.dot(v) * (1 - c) + v.cross(axis) * s;
}

//Returns the position of the particle, after rotating the camera
Vec3f adjParticlePos(Vec3f pos) {
	return rotate(pos, Vec3f(1, 0, 0), -30);
}

//Returns whether particle1 is in back of particle2
bool compareParticles(Particle* particle1, Particle* particle2) {
	return adjParticlePos(particle1->pos)[2] <
		adjParticlePos(particle2->pos)[2];
}


		
ParticleEngine::ParticleEngine(GLuint textureId1) {
	textureId = textureId1;
	timeUntilNextStep = 0;
	angle = 0;
	for(int i = 0; i < NUM_PARTICLES; i++) {
		createParticle(particles + i);
	}
	for(int i = 0; i < 5 / STEP_TIME; i++) {
		step();
	}
}
		
//Advances the particle fountain by the specified amount of time.
void ParticleEngine::advance(float dt) {
	while (dt > 0) {
		if (timeUntilNextStep < dt) {
			dt -= timeUntilNextStep;
			step();
			timeUntilNextStep = STEP_TIME;
		}
		else {
			timeUntilNextStep -= dt;
			dt = 0;
		}
	}
}
		
//Draws the particle fountain.
void ParticleEngine::draw() {
	vector<Particle*> ps;
	for(int i = 0; i < NUM_PARTICLES; i++) {
		ps.push_back(particles + i);
	}
	sort(ps.begin(), ps.end(), compareParticles);
			
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			
	glBegin(GL_QUADS);
	for(unsigned int i = 0; i < ps.size(); i++) {
		Particle* p = ps[i];
		glColor4f(p->color[0], p->color[1], p->color[2],
					(1 - p->timeAlive / p->lifespan));
		float size = p->size / 2;
				
		Vec3f pos = adjParticlePos(p->pos);
				
		glTexCoord2f(0, 0);
		glVertex3f(pos[0] - size, pos[1] - size, pos[2]);
		glTexCoord2f(0, 1);
		glVertex3f(pos[0] - size, pos[1] + size, pos[2]);
		glTexCoord2f(1, 1);
		glVertex3f(pos[0] + size, pos[1] + size, pos[2]);
		glTexCoord2f(1, 0);
		glVertex3f(pos[0] + size, pos[1] - size, pos[2]);
	}
	glEnd();
}

//Returns an array indicating pixel data for an RGBA image that is the same as
//image, but with an alpha channel indicated by the grayscale image alphaChannel
char* addAlphaChannel(Image* image, Image* alphaChannel) {
	char* pixels = new char[image->width * image->height * 4];
	for(int y = 0; y < image->height; y++) {
		for(int x = 0; x < image->width; x++) {
			for(int j = 0; j < 3; j++) {
				pixels[4 * (y * image->width + x) + j] =
					image->pixels[3 * (y * image->width + x) + j];
			}
			pixels[4 * (y * image->width + x) + 3] =
				alphaChannel->pixels[3 * (y * image->width + x)];
		}
	}
	
	return pixels;
}

//Makes the image into a texture, using the specified grayscale image as an
//alpha channel and returns the id of the texture
GLuint loadAlphaTexture(Image* image, Image* alphaChannel) {
	char* pixels = addAlphaChannel(image, alphaChannel);
	
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexImage2D(GL_TEXTURE_2D,
				 0,
				 GL_RGBA,
				 image->width, image->height,
				 0,
				 GL_RGBA,
				 GL_UNSIGNED_BYTE,
				 pixels);
	
	delete pixels;
	return textureId;
}

// draw the particles system
void drawParticles(ParticleEngine* input, double tx, double ty, double tz, double sx, double sy, double sz, R3Point *shipPos, double cullFront, double cullBack) {

	double diff = ty - shipPos->Y();
	if ((diff < cullFront) && (diff > -cullBack)) {
		glPushMatrix();
		glTranslatef(tx, ty, tz);
		glScalef(sx, sy, sz);
		glRotatef(90,1,0,0);
		input->draw();
		glPopMatrix();
	}
}
