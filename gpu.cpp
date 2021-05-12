/*!
 * @file
 * @brief This file contains implementation of gpu
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */

#include <student/gpu.hpp>

#include <stdbool.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

// Minimum ze dvou zadanych hodnot
#define MIN(a, b)       (((a) < (b)) ? (a) : (b))

// Maximum ze dvou zadanych hodnot
#define MAX(a, b)       (((a) > (b)) ? (a) : (b))

struct Triangle{
  OutVertex points[3];
};

/*
struct Triangle {
  OutVertex point[3];
};

InVertex computeVertexID(GPUContext &ctx, InVertex inVertex, uint32_t i) {
    // kontrola, zda je zapnuté indexování
    if(ctx.vao.indexBuffer){  // vao.indexBuffer =! nullptr -> indexování je zapnuto
      // 32-bitové indexování
      if(ctx.vao.indexType == IndexType::UINT32)
        inVertex.gl_VertexID = ((uint32_t*)ctx.vao.indexBuffer)[i];
      // 16-bitové indexování
      else if(ctx.vao.indexType == IndexType::UINT16)
        inVertex.gl_VertexID = ((uint16_t*)ctx.vao.indexBuffer)[i];
      // 8-bitové indexování
      else if(ctx.vao.indexType == IndexType::UINT8 )
        inVertex.gl_VertexID = ((uint8_t *)ctx.vao.indexBuffer)[i];

    }else{  // vao.indexBuffer == nullptr -> indexování není zapnuto
      inVertex.gl_VertexID = i;
    }
    return inVertex;
}
*/
/*
InVertex readAtributes(GPUContext &ctx, InVertex inVertex) {
    for(uint32_t j = 0; j < maxAttributes; j++){
      auto const &currPos = ctx.vao.vertexAttrib[j];
      // pokud nejsou na dané pozici žádná data ke čtení, cyklus pokračuje dál
      if(!currPos.bufferData)
        continue;
      // pokud jsou na dané pozici nějaká data, zjistíme jaká a přečteme je z adresy buf_ptr + offset + stride*gl_VertexID
      if(currPos.type == AttributeType::FLOAT)
        inVertex.attributes[j].v1 = *(float*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
      if(currPos.type == AttributeType::VEC2 )
        inVertex.attributes[j].v2 = *(glm::vec2*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
      if(currPos.type == AttributeType::VEC3 )
        inVertex.attributes[j].v3 = *(glm::vec3*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
      if(currPos.type == AttributeType::VEC4 )
        inVertex.attributes[j].v4 = *(glm::vec4*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
    }
    return inVertex;
}
*/
/*
void runVertexAssembly(GPUContext &ctx, InVertex *inVertex, uint32_t v){
  computeVertexID(ctx, inVertex, v);
  readAtributes(ctx, inVertex);
}
*/


/*
void runPrimitiveAssembly(Triangle *triangle, GPUContext &ctx, uint32_t t) {
  for(uint32_t v = 0; v < 3; v++){
    InVertex inVertex;
    // triangle->points[0].gl_Position.x = 2;
    // runVertexAssembly(ctx, inVertex, t + v);

      inVertex = computeVertexID(ctx, inVertex, v);
      inVertex = readAtributes(ctx, inVertex);

        printf("before: x=%g, y=%g\n", triangle->point[0].gl_Position.x, triangle->point[0].gl_Position.y);

    ctx.prg.vertexShader(triangle->point[v], inVertex, ctx.prg.uniforms);
     printf("after: x=%g, y=%g\n", triangle->point[1].gl_Position.x, triangle->point[1].gl_Position.y);

 
  }
}
*/
/*
void createFragment(InFragment &inFragment, Triangle &triangle, glm::vec3 barycentric, glm::vec2 pixelCoord, GPUContext &ctx) {

}
*/
/*
float TriangleArea(float x0, float y0, float x1, float y1, float x2, float y2)
{
    float area_triangle;
    float a, b, c, s;

    a=std::sqrt(ABS((x0-x1)*(x0-x1)-(y0-y1)*(y0-y1)));
    b=std::sqrt(ABS((x1-x2)*(x1-x2)-(y1-y2)*(y1-y2)));
    c=std::sqrt(ABS((x0-x2)*(x0-x2)-(y0-y2)*(y0-y2)));

    s=(a+b+c)/2;

    area_triangle=std::sqrt(ABS((s*(s-a)*(s-b)*(s-c))));
    
    //printf("area: %g\n", area_triangle);

    return area_triangle;

}
*/

void runPerspectiveDivision(Triangle *triangle) {
  for(unsigned i = 0; i < 3; i++) {
    triangle->points[i].gl_Position.x = triangle->points[i].gl_Position.x / triangle->points[i].gl_Position.w;
    triangle->points[i].gl_Position.y = triangle->points[i].gl_Position.y / triangle->points[i].gl_Position.w;
    triangle->points[i].gl_Position.z = triangle->points[i].gl_Position.z / triangle->points[i].gl_Position.w;
  }
}

void runViewportTransformation(Triangle *triangle, GPUContext &ctx) {
  for (unsigned i = 0; i < 3; i++) {
    triangle->points[i].gl_Position.x = (triangle->points[i].gl_Position.x * 0.5 + 0.5) * (ctx.frame.width);
    triangle->points[i].gl_Position.y = (triangle->points[i].gl_Position.y * 0.5 + 0.5) * (ctx.frame.height);
    // z a w zustava nezmenene po perspective division
  }
}

InVertex computeVertexID(InVertex inVertex, GPUContext &ctx, int i){
  if(ctx.vao.indexBuffer != nullptr){
    if(ctx.vao.indexType == IndexType::UINT8){
      uint8_t *ind = (uint8_t*)ctx.vao.indexBuffer;
      inVertex.gl_VertexID = ind[i];
    }else if(ctx.vao.indexType == IndexType::UINT16){
      uint16_t *ind = (uint16_t*)ctx.vao.indexBuffer;
      inVertex.gl_VertexID = ind[i];
    }else{
      uint32_t *ind = (uint32_t*)ctx.vao.indexBuffer;
      inVertex.gl_VertexID = ind[i];
    }
  }else{
    inVertex.gl_VertexID = i;
  }
  return inVertex;
}

InVertex readAttributes(InVertex inVertex, GPUContext &ctx){
  for(int j = 0; j < maxAttributes; j++){
    if(ctx.vao.vertexAttrib[j].type != AttributeType::EMPTY){
      //Convert buffer to char*, so we can do pointer arithmetic by bytes
      char *buffer = (char*)ctx.vao.vertexAttrib[j].bufferData;
      int offset = ctx.vao.vertexAttrib[j].offset 
        + ctx.vao.vertexAttrib[j].stride * inVertex.gl_VertexID;

      if(ctx.vao.vertexAttrib[j].type == AttributeType::FLOAT){
        inVertex.attributes[j].v1 = *(float*)(buffer + offset);
      }
      else if(ctx.vao.vertexAttrib[j].type == AttributeType::VEC2){
        inVertex.attributes[j].v2 = *(glm::vec2*)(buffer + offset);
      }
      else if(ctx.vao.vertexAttrib[j].type == AttributeType::VEC3){
        inVertex.attributes[j].v3 = *(glm::vec3*)(buffer + offset);
      }
      else if(ctx.vao.vertexAttrib[j].type == AttributeType::VEC4){
        inVertex.attributes[j].v4 = *(glm::vec4*)(buffer + offset);
      }
    }
  }
  return inVertex;
}

void loadTriangle(GPUContext &ctx, uint32_t nofVertices, Triangle *triangle, int t){
  for(int i = t; i < t + 3; i++){//smyčka přes vrcholy
    InVertex inVertex; // vrchol, co vstupuje do vertex shader

    inVertex = computeVertexID(inVertex, ctx, i);
    inVertex = readAttributes(inVertex, ctx);

    ctx.prg.vertexShader(triangle->points[i - t],inVertex,ctx.prg.uniforms);
  }
}

float dist(float p1[], float p2[]){
  return (float)sqrt(pow((double)(p2[0] - p1[0]), 2) + pow((double)(p2[1] - p1[1]), 2) * 1.0f);
}

float area(float a, float b, float c){
  float s = (a+b+c)/2;
  s = (s*(s-a)*(s-b)*(s-c));
  return (float)sqrt((double)((s < 0) ? -s : s));
}

glm::vec3 calcBary2(float p[], float a[], float b[], float c[]){
  float abc = area(dist(a, b), dist(b, c), dist(c, a));
  float bcp = area(dist(b, c), dist(c, p), dist(p, b));
  float acp = area(dist(a, c), dist(c, p), dist(p, a));
  float abp = area(dist(a, b), dist(b, p), dist(p, a));
  return glm::vec3(bcp/abc, acp/abc, abp/abc);
}

glm::vec3 getBary(Triangle *triangle, float x, float y){
  float p[4][2];
  for(int i = 0; i < 2; i++){
    p[0][i] = (float)triangle->points[0].gl_Position[i];
    p[1][i] = (float)triangle->points[1].gl_Position[i];
    p[2][i] = (float)triangle->points[2].gl_Position[i];
  }
  p[3][0] = x;
  p[3][1] = y;
  return glm::vec3(calcBary2(p[3], p[0], p[1], p[2]));
}

void rasterize(GPUContext &ctx, Triangle *triangle) {
  // spočítat hranice trojúhelníku
  int minX = MIN(MIN(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int minY = MIN(MIN(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);
  int maxX = MAX(MAX(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int maxY = MAX(MAX(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);

	minX = MAX(0, minX);
	minY = MAX(0, minY);
	maxX = MIN(ctx.frame.width - 1, maxX);
	maxY = MIN(ctx.frame.height - 1, maxY);

  // cyklus projde celý frame buffer
	for (int y = minY; y <= maxY; y++) {
		bool even = (y - minY) % 2 == 0;

		int startX = even ? minX : maxX;
		int endX = even ? maxX : minX;
		int stepX = even ? 1 : -1;

		for (int x = startX; x != endX; x += stepX) {
        if (even && x > endX || !even && x < endX) {
          break;
        }
        if (x != ctx.frame.width && y != ctx.frame.height) {
          glm::vec3 bary = getBary(triangle, x+0.5f, y+0.5f);
          float lambda_sum = bary[0] + bary[1] + bary[2];
          if (lambda_sum < 1.001f && lambda_sum > 0.999f) {
            InFragment inFragment;
            inFragment.gl_FragCoord = glm::vec4((double)x + 0.5f, (double)y + 0.5f, 1.0f, 1.0f);
            OutFragment outFragment;
            ctx.prg.fragmentShader(outFragment, inFragment, ctx.prg.uniforms);
          }
        }
      }
  }
}

//! [drawTrianglesImpl]
void drawTrianglesImpl(GPUContext &ctx,uint32_t nofVertices){
  /// \todo Tato funkce vykreslí trojúhelníky podle daného nastavení.<br>
  /// ctx obsahuje aktuální stav grafické karty.
  /// Parametr "nofVertices" obsahuje počet vrcholů, který by se měl vykreslit (3 pro jeden trojúhelník).<br>
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace.

Triangle triangle;
for(int i = 0; i < nofVertices; i += 3){
  loadTriangle(ctx, nofVertices, &triangle, i);
    runPerspectiveDivision(&triangle);
    runViewportTransformation(&triangle, ctx);
    rasterize(ctx, &triangle);
  }
}

//! [drawTrianglesImpl]

/**
 * @brief This function reads color from texture.
 *
 * @param texture texture
 * @param uv uv coordinates
 *
 * @return color 4 floats
 */
glm::vec4 read_texture(Texture const&texture,glm::vec2 uv){
  if(!texture.data)return glm::vec4(0.f);
  auto uv1 = glm::fract(uv);
  auto uv2 = uv1*glm::vec2(texture.width-1,texture.height-1)+0.5f;
  auto pix = glm::uvec2(uv2);
  //auto t   = glm::fract(uv2);
  glm::vec4 color = glm::vec4(0.f,0.f,0.f,1.f);
  for(uint32_t c=0;c<texture.channels;++c)
    color[c] = texture.data[(pix.y*texture.width+pix.x)*texture.channels+c]/255.f;
  return color;
}

/**
 * @brief This function clears framebuffer.
 *
 * @param ctx GPUContext
 * @param r red channel
 * @param g green channel
 * @param b blue channel
 * @param a alpha channel
 */
void clear(GPUContext&ctx,float r,float g,float b,float a){
  auto&frame = ctx.frame;
  auto const nofPixels = frame.width * frame.height;
  for(size_t i=0;i<nofPixels;++i){
    frame.depth[i] = 10e10f;
    frame.color[i*4+0] = static_cast<uint8_t>(glm::min(r*255.f,255.f));
    frame.color[i*4+1] = static_cast<uint8_t>(glm::min(g*255.f,255.f));
    frame.color[i*4+2] = static_cast<uint8_t>(glm::min(b*255.f,255.f));
    frame.color[i*4+3] = static_cast<uint8_t>(glm::min(a*255.f,255.f));
  }
}

