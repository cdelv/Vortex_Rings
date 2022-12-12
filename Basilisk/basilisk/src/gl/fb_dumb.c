/**
# "Dumb" OpenGL implementation

This provides the minimal set of OpenGL functions necessary to
compile/run the interactive version of Basilisk View. It has no
dependencies other than the standard C library. */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "gl.h"
#include "framebuffer.h"

/**
## "Unused" OpenGL functions */

void glBegin (GLenum mode) {}
void glEnd (void) {}
void glClear (GLbitfield mask) {}
void glClearColor (float red, float green, float blue, float alpha) {}
void glBindTexture (GLenum target, GLuint texture) {}
void glColor3f (GLfloat red, GLfloat green, GLfloat blue) {}
void glColorMaterial (GLenum face, GLenum mode) {}
void glDisable (GLenum cap) {}
void glEnable (GLenum cap) {}
void glFinish (void) {}
void glGetDoublev (GLenum pname, GLdouble * params) {}
void glGetIntegerv (GLenum pname, GLint * params) {}
void glHint (GLenum target, GLenum mode) {}
void glLightfv (GLenum light, GLenum pname, const GLfloat *params) {}
void glLightModeli (GLenum pname, GLint param) {}
void glLineWidth (GLfloat width) {}
void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz) {}
void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top,
	      GLdouble nearVal, GLdouble farVal) {}
void glShadeModel (GLenum mode) {}
void glTexCoord1d (GLdouble s) {}
void glTexCoord2f (GLfloat s, GLfloat t) {}
void glTexImage1D (GLenum target, GLint level, GLint internalFormat,
		   GLsizei width, GLint border, GLenum format, GLenum type,
		   const void * data) {}
void glTexParameteri (GLenum target, GLenum pname, GLint param) {}
void gluPerspective (GLdouble fovy, GLdouble aspect, GLdouble zNear,
		     GLdouble zFar) {}
void glVertex3d (GLdouble x, GLdouble y, GLdouble z) {}
void glVertex3f (GLfloat x, GLfloat y, GLfloat z) {}

/**
## Matrix transformations */

static GLfloat modelview[16] = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
}, * current_modelview = modelview;

static GLfloat current_other[16] = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
};

static GLfloat * current = current_other;

GLenum glGetError (void) { return GL_NO_ERROR; }

void glGetFloatv (GLenum pname, GLfloat * params)
{
  int i;
  if (pname == GL_MODELVIEW_MATRIX)
    for (i = 0; i < 16; i++)
      params[i] = current_modelview[i];
}

void glMultMatrixf (const GLfloat * m)
{
  GLfloat r[16];
  int i, j, k;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      GLfloat a = 0.;
      for (k = 0; k < 4; k++)
	a += current[4*k + i]*m[4*j + k];
      r[4*j + i] = a;
    }
  for (i = 0; i < 16; i++)
    current[i] = r[i];
}

void glLoadIdentity (void)
{
  static GLfloat identity[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  };
  int i;
  for (i = 0; i < 16; i++)
    current[i] = identity[i];
}

void glScalef (GLfloat x, GLfloat y, GLfloat z)
{
  glMultMatrixf ((GLfloat []){
      x, 0, 0, 0,
      0, y, 0, 0,
      0, 0, z, 0,
      0, 0, 0, 1
  });
}

void glTranslatef (GLfloat x, GLfloat y, GLfloat z)
{
  glMultMatrixf ((GLfloat []){
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      x, y, z, 1
  });
}

void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
  angle *= M_PI/180.;
  float c = cos(angle), s = sin(angle), n = sqrt(x*x + y*y + z*z);
  if (n > 0.)
    x /= n, y /= n, z /= n;
  glMultMatrixf ((GLfloat []){
      x*x*(1. - c) + c,    y*x*(1. - c) + z*s,  x*z*(1. - c) - y*s,   0,
      x*y*(1. - c) - z*s,  y*y*(1. - c) + c,    y*z*(1. - c) + x*s,   0,
      x*z*(1. - c) + y*s,  y*z*(1. - c) - x*s,  z*z*(1. - c) + c,     0,
      0,                   0,                   0,                    1
  });  
}

void glMatrixMode (GLenum mode)
{
  if (mode == GL_MODELVIEW)
    current = current_modelview;
  else
    current = current_other;
}

static void * stack = NULL;
static int len = 0;

#define SIZE (16*sizeof (GLfloat))

void glPopMatrix (void)
{
  if (current == current_modelview) {
    len -= SIZE;
    assert (len >= 0);
    if (len == 0) {
      free (stack);
      stack = NULL;
      current = current_modelview = modelview;
    }
    else {
      stack = realloc (stack, len);
      current = current_modelview = (GLfloat *)(stack + len - SIZE);
    }
  } 
}

void glPushMatrix (void)
{
  if (current == current_modelview) {
    stack = realloc (stack, len + SIZE);
    current_modelview = (GLfloat *)(stack + len);
    memcpy (current_modelview, current, SIZE);
    current = current_modelview;
    len += SIZE;
  }  
}

/**
## Framebuffer */
    
struct _framebuffer {
};

void framebuffer_destroy (framebuffer * p) {}

framebuffer * framebuffer_new (unsigned width, unsigned height)
{
  return NULL;
}

unsigned char * framebuffer_image (framebuffer * p) {
  return NULL;
}

fbdepth_t * framebuffer_depth (framebuffer * p) {
  return NULL;
}
