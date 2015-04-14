#
# Makefile for producing Assignment #0 output
#

EXE=src/particleview

IMGS=	clean	\
output	\
art/art	\
output/lifetimes	\
output/sinks	\
output/euler	\
output/sources	\
output/boids	\
output/interactions	\


all: $(IMGS)
output:
	mkdir $@

output/interactions: input/
	mkdir $@
	$(EXE) $^opposites.scn -video_prefix $@/opposites -recordandquit 200 -slow
	convert -delay 4 -loop 0 $@/opposites*.jpg output/mutualparticleattraction.gif
	$(EXE) $^rope.scn -video_prefix $@/rope -recordandquit 500 -slow
	convert -delay 4 -loop 0 $@/rope*.jpg output/rope.gif
	$(EXE) $^cloth.scn -video_prefix $@/cloth -recordandquit 300 -slow
	convert -delay 4 -loop 0 $@/cloth*.jpg output/cloth.gif
	
output/lifetimes: input/differentlifetimes.scn
	mkdir $@
	$(EXE) $^ -video_prefix $@/lifetimes -recordandquit 200
	convert -delay 4 -loop 0 $@/lifetimes*.jpg output/differentlifetimes.gif

output/sinks: input/
	mkdir $@
	$(EXE) $^spheresink.scn -video_prefix $@/sphere -recordandquit 500
	convert -delay 4 -loop 0 $@/sphere*.jpg output/spheresink.gif
	$(EXE) $^meshsink.scn -video_prefix $@/mesh -recordandquit 500
	convert -delay 4 -loop 0 $@/mesh*.jpg output/meshsink.gif
	$(EXE) $^boxsink.scn -video_prefix $@/box -recordandquit 500
	convert -delay 4 -loop 0 $@/box*.jpg output/boxsink.gif

output/euler: input/
	mkdir $@
	$(EXE) $^oneparticle.scn -video_prefix $@/oneparticle -recordandquit 500
	convert -delay 4 -loop 0 $@/oneparticle*.jpg output/particlemotion.gif
	$(EXE) $^cannonball.scn -video_prefix $@/cannonball -recordandquit 100
	convert -delay 4 -loop 0 $@/cannonball*.jpg output/effectofgravity.gif
	$(EXE) $^dragrace.scn -video_prefix $@/drag -recordandquit 500
	convert -delay 4 -loop 0 $@/drag*.jpg output/effectofdrag.gif
	$(EXE) $^collisiondetection.scn -video_prefix $@/collision -recordandquit 500
	convert -delay 4 -loop 0 $@/collision*.jpg output/collisiondetection.gif

output/sources: input/sources.scn
	mkdir $@
	$(EXE) $^ -video_prefix $@/sources -recordandquit 500
	convert -delay 4 -loop 0 $@/sources*.jpg output/sources.gif

output/boids: input/flocking.scn
	mkdir $@
	$(EXE) $^ -video_prefix $@/boids -recordandquit 500 -boid
	convert -delay 4 -loop 0 $@/boids*.jpg output/flocking.gif

art/art: input/art.scn
	mkdir $@
	$(EXE) $^ -video_prefix $@/art -recordandquit 100
	convert -delay 4 -loop 0 $@/art*.jpg art/LUMOS.gif
	
clean:
	rm -f -r output



