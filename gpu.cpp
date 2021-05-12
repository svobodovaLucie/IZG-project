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

InVertex computeVertexID(InVertex inVertex, GPUContext &ctx, uint32_t i) {
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

InVertex readAttributes(InVertex inVertex, GPUContext &ctx) {
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

glm::vec3 get_bary(Triangle *triangle, float x, float y){
  float p[4][2];
  for(int i = 0; i < 2; i++){
    p[0][i] = (float)triangle->points[0].gl_Position[i];
    p[1][i] = (float)triangle->points[1].gl_Position[i];
    p[2][i] = (float)triangle->points[2].gl_Position[i];
  }
  p[3][0] = x;
  p[3][1] = y;

  float abc = area(dist(p[0], p[1]), dist(p[1], p[2]), dist(p[2], p[0]));
  return glm::vec3(area(dist(p[1], p[2]), dist(p[2], p[3]), dist(p[3], p[1]))/abc, area(dist(p[0], p[2]), dist(p[2], p[3]), dist(p[3], p[0]))/abc, area(dist(p[0], p[1]), dist(p[1], p[3]), dist(p[3], p[0]))/abc);
}

void rasterize(GPUContext &ctx, Triangle *triangle) {
  // spočítat hranice trojúhelníku
  int minX = MIN(MIN(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int minY = MIN(MIN(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);
  int maxX = MAX(MAX(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int maxY = MAX(MAX(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);

	minX = MAX(0, minX);
	minY = MAX(0, minY);
	maxX = MIN(ctx.frame.width, maxX);
	maxY = MIN(ctx.frame.height, maxY);

  // cyklus projde celý frame buffer
	for (int y = minY; y <= maxY; y++) {
		bool even = (y - minY) % 2 == 0;

		int startX = even ? minX : maxX;
		int endX = even ? maxX : minX;
		int stepX = even ? 1 : -1;

    // cyklus projde celou
		for (int x = startX; x != endX; x += stepX) {
      if (even && x > endX || !even && x < endX) {
        break;
      }
      if (x != ctx.frame.width && y != ctx.frame.height) {
        glm::vec3 bary = get_bary(triangle, x+0.5f, y+0.5f);
        float lambda_sum = bary[0] + bary[1] + bary[2];
        if (lambda_sum < 1.001f && lambda_sum > 0.999f) {
          float depth_before = triangle->points[0].gl_Position.z * bary[0] + triangle->points[1].gl_Position.z * bary[1] + triangle->points[2].gl_Position.z * bary[2];
          InFragment inFragment;
          // interpolace
          float s = bary[0]/triangle->points[0].gl_Position.w + bary[1]/triangle->points[1].gl_Position.w + bary[2]/triangle->points[2].gl_Position.w;
          bary[0] = bary[0] / (triangle->points[0].gl_Position.w*s);
          bary[1] = bary[1] / (triangle->points[1].gl_Position.w*s);
          bary[2] = bary[2] / (triangle->points[2].gl_Position.w*s);
          for(int i = 0; i < maxAttributes; i++){
            if(ctx.prg.vs2fs[i] != AttributeType::EMPTY)
              inFragment.attributes[i].v4 = triangle->points[0].attributes[i].v4 * bary[0] + triangle->points[1].attributes[i].v4 * bary[1] + triangle->points[2].attributes[i].v4 * bary[2];
          }
          // fragment
          inFragment.gl_FragCoord = glm::vec4((double)x + 0.5f, (double)y + 0.5f, depth_before, 1.0f);
          OutFragment outFragment;
          ctx.prg.fragmentShader(outFragment, inFragment, ctx.prg.uniforms);
        
          // zapis do frame bufferu 
          if(ctx.frame.depth[x + y * ctx.frame.width] > inFragment.gl_FragCoord.z){
            *(x * 4 + ctx.frame.color+ 4 * y * ctx.frame.width) = ((float)*(x * 4 + ctx.frame.color + 4 * y * ctx.frame.width)/255.f * (1.0f-outFragment.gl_FragColor[3])+outFragment.gl_FragColor[0]*outFragment.gl_FragColor[3])*255.f;
            *(x * 4 + 1 + ctx.frame.color + 4 * y * ctx.frame.width) = ((float)*(x * 4 + 1 + ctx.frame.color + 4 * y * ctx.frame.width)/255.f * (1.0f-outFragment.gl_FragColor[3])+outFragment.gl_FragColor[1]*outFragment.gl_FragColor[3])*255.f;
            *(x * 4 + 2 + ctx.frame.color + 4 * y * ctx.frame.width) = ((float)*(x * 4 + 2 + ctx.frame.color + 4 * y * ctx.frame.width)/255.f * (1.0f-outFragment.gl_FragColor[3])+outFragment.gl_FragColor[2]*outFragment.gl_FragColor[3])*255.f;
            if (outFragment.gl_FragColor[3]>0.5f) {
                *(x + ctx.frame.depth + y * ctx.frame.width) = inFragment.gl_FragCoord.z;
            }
          }
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

